#!/usr/bin/env python
"""Combine STRING disease-module matching with TDP phenotype specificity.

No OpenTargets or expression-matrix evidence is used.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--string-scores",
        type=Path,
        default=root / "outputs" / "string_module_disease_matching" / "string_module_disease_scores.csv",
    )
    parser.add_argument(
        "--tdp-scores",
        type=Path,
        default=root / "outputs" / "tdp_phenotype_matching" / "tdp_match_scores.csv",
    )
    parser.add_argument("--out-dir", type=Path, default=root / "outputs" / "string_tdp_final")
    parser.add_argument("--string-weight", type=float, default=0.65)
    parser.add_argument("--tdp-weight", type=float, default=0.35)
    args = parser.parse_args()

    string_rows = {row["disease"]: row for row in read_csv(args.string_scores)}
    tdp_rows = {row["disease"]: row for row in read_csv(args.tdp_scores)}
    diseases = sorted(set(string_rows) | set(tdp_rows))
    rows: list[dict[str, Any]] = []
    for disease in diseases:
        string_row = string_rows.get(disease, {})
        tdp_row = tdp_rows.get(disease, {})
        string_score = safe_float(string_row.get("final_string_module_disease_score"))
        tdp_score = safe_float(tdp_row.get("tdp_phenotype_match_score"))
        final = args.string_weight * string_score + args.tdp_weight * tdp_score
        rows.append(
            {
                "disease": disease,
                "final_string_tdp_score": round(final, 4),
                "string_module_disease_score": round(string_score, 4),
                "tdp_phenotype_match_score": round(tdp_score, 4),
                "string_specificity_overlap_score": string_row.get("specificity_weighted_overlap_score", ""),
                "string_pathway_axis_score": string_row.get("string_pathway_axis_score", ""),
                "tdp_top_matched_axes": tdp_row.get("top_matched_axes", ""),
                "string_overlap_genes": string_row.get("overlap_genes", ""),
            }
        )

    rows.sort(key=lambda row: row["final_string_tdp_score"], reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank

    fields = [
        "rank",
        "disease",
        "final_string_tdp_score",
        "string_module_disease_score",
        "tdp_phenotype_match_score",
        "string_specificity_overlap_score",
        "string_pathway_axis_score",
        "tdp_top_matched_axes",
        "string_overlap_genes",
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "final_string_tdp_scores.csv", rows, fields)

    md = [
        "# Final STRING + TDP Autoimmune Indication Scores",
        "",
        "No OpenTargets or expression-matrix evidence is used.",
        "",
        "```text",
        f"Final = {args.string_weight:.0%} STRING disease-module score + {args.tdp_weight:.0%} TDP phenotype specificity score",
        "```",
        "",
        "| Rank | Disease | Final | STRING module | TDP phenotype | Key TDP axes |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for row in rows:
        axes = str(row["tdp_top_matched_axes"]).replace(" | ", "<br>")
        md.append(
            f"| {row['rank']} | {row['disease']} | {row['final_string_tdp_score']:.2f} | "
            f"{row['string_module_disease_score']:.2f} | {row['tdp_phenotype_match_score']:.2f} | "
            f"{axes} |"
        )
    (args.out_dir / "final_string_tdp_scores.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print({"ok": True, "rows": len(rows), "top_disease": rows[0]["disease"], "out_dir": str(args.out_dir)})


if __name__ == "__main__":
    main()
