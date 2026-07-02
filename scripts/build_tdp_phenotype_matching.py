#!/usr/bin/env python
"""Target degradation phenotype matching without expression-matrix downloads.

The workflow is a deterministic prototype for the user's TDP concept:

1. Convert measured target degradation weights into a phenotype-axis profile.
2. Match that TDP profile against disease phenotype profiles.
3. Export ranked indication scores and axis-level contributions.

No Open Targets or expression-matrix evidence is used.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


CELL_TYPE_AXES = {
    "b_cell_receptor_autoantibody",
    "mast_basophil_fceri",
    "t_cell_receptor_signaling",
    "myeloid_tlr_tnf_il6",
    "epithelial_barrier_signaling",
    "type_i_interferon",
    "th2_cytokine",
    "il23_il17_tissue_inflammation",
}


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def load_target_weights(path: Path) -> dict[str, float]:
    payload = read_json(path)
    weights: dict[str, float] = {}
    for row in payload.get("targets", []):
        symbol = str(row.get("symbol", "")).upper().strip()
        if symbol and not symbol.startswith("REPLACE_WITH"):
            weights[symbol] = safe_float(row.get("weight"))
    if not weights:
        raise ValueError(
            "No usable target weights found. Provide a local target profile with real human gene symbols."
        )
    return weights


def build_tdp_axis_profile(
    target_weights: dict[str, float],
    target_axis_effects: dict[str, dict[str, float]],
    axes: list[str],
) -> tuple[dict[str, float], list[dict[str, Any]]]:
    raw = {axis: 0.0 for axis in axes}
    contribution_rows: list[dict[str, Any]] = []

    for target, degradation_weight in target_weights.items():
        for axis, target_axis_effect in target_axis_effects.get(target, {}).items():
            contribution = degradation_weight * safe_float(target_axis_effect)
            raw[axis] = raw.get(axis, 0.0) + contribution
            contribution_rows.append(
                {
                    "target": target,
                    "target_degradation_weight": round(degradation_weight, 4),
                    "axis": axis,
                    "target_axis_effect_prior": round(safe_float(target_axis_effect), 4),
                    "raw_axis_contribution": round(contribution, 6),
                }
            )

    max_raw = max(raw.values()) if raw else 0.0
    if max_raw <= 0:
        normalized = {axis: 0.0 for axis in axes}
    else:
        normalized = {axis: value / max_raw * 100.0 for axis, value in raw.items()}

    return normalized, contribution_rows


def weighted_axis_coverage(
    disease_profile: dict[str, float],
    tdp_axis_score: dict[str, float],
    axes: set[str] | None = None,
) -> float:
    total_weight = 0.0
    weighted_sum = 0.0
    for axis, need in disease_profile.items():
        if axes is not None and axis not in axes:
            continue
        weight = safe_float(need)
        total_weight += weight
        weighted_sum += weight * safe_float(tdp_axis_score.get(axis))
    return weighted_sum / total_weight if total_weight > 0 else 0.0


def rank_diseases(
    tdp_axis_score: dict[str, float],
    disease_profiles: dict[str, dict[str, float]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    disease_rows: list[dict[str, Any]] = []
    axis_rows: list[dict[str, Any]] = []

    for disease, profile in disease_profiles.items():
        axis_coverage = weighted_axis_coverage(profile, tdp_axis_score)
        cell_type_coverage = weighted_axis_coverage(profile, tdp_axis_score, CELL_TYPE_AXES)
        final = 0.80 * axis_coverage + 0.20 * cell_type_coverage

        axis_contribs: list[tuple[str, float]] = []
        total_need = sum(safe_float(v) for v in profile.values())
        for axis, need in profile.items():
            need_value = safe_float(need)
            tdp_value = safe_float(tdp_axis_score.get(axis))
            contribution = need_value * tdp_value / total_need if total_need > 0 else 0.0
            axis_contribs.append((axis, contribution))
            axis_rows.append(
                {
                    "disease": disease,
                    "axis": axis,
                    "disease_axis_need_0_1": round(need_value, 4),
                    "tdp_axis_activity_0_100": round(tdp_value, 4),
                    "axis_match_contribution": round(contribution, 4),
                }
            )

        axis_contribs.sort(key=lambda item: item[1], reverse=True)
        top_axes = " | ".join(f"{axis}:{score:.2f}" for axis, score in axis_contribs[:5])

        disease_rows.append(
            {
                "disease": disease,
                "tdp_phenotype_match_score": round(final, 4),
                "axis_coverage_score": round(axis_coverage, 4),
                "cell_type_specificity_score": round(cell_type_coverage, 4),
                "top_matched_axes": top_axes,
                "disease_axis_count": len(profile),
                "tdp_data_completeness_note": "Measured depth for targets; duration/resynthesis/scaffold removal/omics shifts are proxy-only until assay data are supplied.",
            }
        )

    disease_rows.sort(key=lambda row: row["tdp_phenotype_match_score"], reverse=True)
    for rank, row in enumerate(disease_rows, start=1):
        row["rank"] = rank

    return disease_rows, axis_rows


def write_markdown(path: Path, disease_rows: list[dict[str, Any]], tdp_axis_rows: list[dict[str, Any]]) -> None:
    lines: list[str] = []
    lines.append("# Target degradation phenotype matching")
    lines.append("")
    lines.append("No Open Targets or expression-matrix evidence was used.")
    lines.append("")
    lines.append("## Formula")
    lines.append("")
    lines.append("```text")
    lines.append("TDP axis activity = sum(target degradation weight * target-axis effect prior)")
    lines.append("axis coverage = sum(disease axis need * TDP axis activity) / sum(disease axis need)")
    lines.append("cell-type specificity = same calculation restricted to cell-type axes")
    lines.append("final TDP match = 80% axis coverage + 20% cell-type specificity")
    lines.append("```")
    lines.append("")
    lines.append("## Ranked indications")
    lines.append("")
    lines.append("| Rank | Disease | TDP match | Axis coverage | Cell-type specificity | Top matched axes |")
    lines.append("|---:|---|---:|---:|---:|---|")
    for row in disease_rows:
        lines.append(
            f"| {row['rank']} | {row['disease']} | {row['tdp_phenotype_match_score']:.2f} | "
            f"{row['axis_coverage_score']:.2f} | {row['cell_type_specificity_score']:.2f} | "
            f"{row['top_matched_axes']} |"
        )
    lines.append("")
    lines.append("## PROTAC TDP axis profile")
    lines.append("")
    lines.append("| Axis | TDP activity |")
    lines.append("|---|---:|")
    for row in tdp_axis_rows:
        lines.append(f"| {row['axis']} | {row['tdp_axis_activity_0_100']:.2f} |")
    lines.append("")
    lines.append("## Interpretation boundary")
    lines.append("")
    lines.append("This is a phenotype-matching model, not a trained efficacy predictor. It becomes stronger when measured phosphoproteome, transcriptome, cytokine, and cell-state shift data are substituted for the current target-mechanism priors.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-profile", type=Path, default=root / "config" / "target_profile.template.json")
    parser.add_argument("--profile-config", type=Path, default=root / "config" / "tdp_disease_profiles.template.json")
    parser.add_argument("--out-dir", type=Path, default=root / "outputs" / "tdp_phenotype_matching")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    target_weights = load_target_weights(args.target_profile)
    config = read_json(args.profile_config)
    axes = list(config["axes"])
    if not config.get("disease_profiles") or any(
        str(name).startswith("REPLACE_WITH") for name in config["disease_profiles"]
    ):
        raise ValueError(
            "No usable disease profiles found. Provide a local profile config with disease phenotype axes."
        )
    tdp_axis_score, target_axis_rows = build_tdp_axis_profile(
        target_weights=target_weights,
        target_axis_effects=config["target_axis_effects"],
        axes=axes,
    )
    disease_rows, disease_axis_rows = rank_diseases(tdp_axis_score, config["disease_profiles"])

    tdp_axis_rows = [
        {"axis": axis, "tdp_axis_activity_0_100": round(score, 4)}
        for axis, score in sorted(tdp_axis_score.items(), key=lambda item: item[1], reverse=True)
    ]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        args.out_dir / "tdp_match_scores.csv",
        disease_rows,
        [
            "rank",
            "disease",
            "tdp_phenotype_match_score",
            "axis_coverage_score",
            "cell_type_specificity_score",
            "top_matched_axes",
            "disease_axis_count",
            "tdp_data_completeness_note",
        ],
    )
    write_csv(args.out_dir / "tdp_axis_profile.csv", tdp_axis_rows, ["axis", "tdp_axis_activity_0_100"])
    write_csv(
        args.out_dir / "target_axis_contributions.csv",
        target_axis_rows,
        [
            "target",
            "target_degradation_weight",
            "axis",
            "target_axis_effect_prior",
            "raw_axis_contribution",
        ],
    )
    write_csv(
        args.out_dir / "disease_axis_contributions.csv",
        disease_axis_rows,
        [
            "disease",
            "axis",
            "disease_axis_need_0_1",
            "tdp_axis_activity_0_100",
            "axis_match_contribution",
        ],
    )
    (args.out_dir / "tdp_match_scores.json").write_text(
        json.dumps(disease_rows, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_markdown(args.out_dir / "tdp_phenotype_matching.md", disease_rows, tdp_axis_rows)

    print(
        json.dumps(
            {
                "ok": True,
                "diseases": len(disease_rows),
                "axes": len(axes),
                "top_disease": disease_rows[0]["disease"] if disease_rows else "",
                "out_dir": str(args.out_dir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
