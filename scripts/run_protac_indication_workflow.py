#!/usr/bin/env python
"""Run the generic PROTAC indication workflow from target/off-target data only.

Input:
  - a target degradation profile with human gene symbols and normalized effect weights

Output:
  - STRING affected protein module
  - Enrichr disease-module provenance
  - STRING/network medicine disease scores
  - automatically inferred TDP phenotype-axis scores
  - final indication ranking
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_string_module_disease_matching import (  # noqa: E402
    ENRICHR_LIBRARIES,
    build_string_module,
    fetch_enrichr_library,
    fetch_string_partners,
    hypergeom_tail,
    load_targets,
    safe_float,
    write_csv,
)


ONCOLOGY_PATTERN = re.compile(
    r"\b(cancer|carcinoma|sarcoma|melanoma|leukemia|lymphoma|myeloma|glioma|blastoma|neoplasm|tumou?r)\b",
    re.IGNORECASE,
)

AXIS_MARKERS = {
    "b_cell_receptor_autoantibody": {
        "CD19",
        "CD79A",
        "CD79B",
        "MS4A1",
        "BLK",
        "BANK1",
        "TNFRSF13B",
        "TNFRSF13C",
        "IGHG1",
        "IGHM",
    },
    "t_cell_receptor_signaling": {
        "CD3D",
        "CD3E",
        "CD3G",
        "CD4",
        "CD8A",
        "CD8B",
        "LCK",
        "ZAP70",
        "LAT",
        "TRAC",
    },
    "mast_basophil_fceri": {
        "FCER1A",
        "FCER1G",
        "KIT",
        "MS4A2",
        "TPSAB1",
        "TPSB2",
        "CPA3",
        "HDC",
        "IL4",
        "IL13",
    },
    "myeloid_tlr_tnf_il6": {
        "TLR2",
        "TLR4",
        "MYD88",
        "NFKB1",
        "RELA",
        "TNF",
        "IL6",
        "IL1B",
        "NLRP3",
        "CXCL8",
    },
    "type_i_interferon": {
        "IFNA1",
        "IFNB1",
        "IRF3",
        "IRF5",
        "IRF7",
        "STAT1",
        "STAT2",
        "ISG15",
        "MX1",
        "OAS1",
    },
    "th2_cytokine": {
        "IL4",
        "IL5",
        "IL13",
        "IL4R",
        "STAT6",
        "GATA3",
        "CCL17",
        "CCL22",
        "TSLP",
        "CRLF2",
    },
    "il23_il17_tissue_inflammation": {
        "IL17A",
        "IL17F",
        "IL23A",
        "IL23R",
        "RORC",
        "STAT3",
        "CCR6",
        "CXCL1",
        "CXCL2",
        "S100A8",
    },
    "epithelial_barrier_signaling": {
        "KRT5",
        "KRT14",
        "KRT17",
        "FLG",
        "LOR",
        "CLDN1",
        "EGF",
        "AREG",
        "EREG",
        "MUC2",
    },
}

DEFAULT_DISEASE_KEYWORDS = [
    "allergy",
    "allergic",
    "arthritis",
    "asthma",
    "atopic",
    "autoimmune",
    "colitis",
    "dermatitis",
    "eczema",
    "immune",
    "inflammation",
    "inflammatory",
    "lupus",
    "psoriasis",
    "sclerosis",
    "spondylitis",
    "urticaria",
]


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def normalize_term(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", " ", value or "")
    return re.sub(r"\s+", " ", text).strip()


def should_keep_disease(term: str, exclude_oncology: bool, disease_keywords: list[str]) -> bool:
    if exclude_oncology and ONCOLOGY_PATTERN.search(term):
        return False
    normalized = normalize_term(term)
    if len(normalized) < 4 or len(normalized.split()) > 12:
        return False
    if not disease_keywords or disease_keywords == ["all"]:
        return True
    normalized_lower = normalized.lower()
    return any(keyword.lower() in normalized_lower for keyword in disease_keywords)


def build_candidate_modules(
    libraries: dict[str, dict[str, set[str]]],
    exclude_oncology: bool,
    disease_keywords: list[str],
    max_candidate_terms: int,
) -> tuple[dict[str, set[str]], list[dict[str, Any]], set[str]]:
    modules: dict[str, set[str]] = {}
    provenance: list[dict[str, Any]] = []
    background: set[str] = set()
    for library, gene_sets in libraries.items():
        for term, genes in gene_sets.items():
            background.update(genes)
            if not should_keep_disease(term, exclude_oncology, disease_keywords):
                continue
            disease = normalize_term(term)
            if not disease:
                continue
            module = modules.setdefault(disease, set())
            before = len(module)
            module.update(genes)
            provenance.append(
                {
                    "disease": disease,
                    "library": library,
                    "matched_term": term,
                    "term_gene_count": len(genes),
                    "new_genes_added": len(module) - before,
                }
            )
            if len(modules) >= max_candidate_terms:
                break
        if len(modules) >= max_candidate_terms:
            break
    return modules, provenance, background


def axis_profile_from_genes(gene_weights: dict[str, float] | set[str]) -> dict[str, float]:
    if isinstance(gene_weights, set):
        weights = {gene: 1.0 for gene in gene_weights}
    else:
        weights = gene_weights

    raw: dict[str, float] = {}
    for axis, markers in AXIS_MARKERS.items():
        raw[axis] = sum(safe_float(weights.get(gene, 0.0)) for gene in markers)

    max_raw = max(raw.values(), default=0.0)
    if max_raw <= 0:
        return {axis: 0.0 for axis in AXIS_MARKERS}
    return {axis: value / max_raw * 100.0 for axis, value in raw.items()}


def cosine_similarity(a: dict[str, float], b: dict[str, float]) -> float:
    axes = sorted(set(a) | set(b))
    dot = sum(safe_float(a.get(axis)) * safe_float(b.get(axis)) for axis in axes)
    a_norm = math.sqrt(sum(safe_float(a.get(axis)) ** 2 for axis in axes))
    b_norm = math.sqrt(sum(safe_float(b.get(axis)) ** 2 for axis in axes))
    if a_norm <= 0 or b_norm <= 0:
        return 0.0
    return dot / (a_norm * b_norm) * 100.0


def score_modules(
    affected_weights: dict[str, float],
    target_weights: dict[str, float],
    disease_modules: dict[str, set[str]],
    background: set[str],
    max_ranked_diseases: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, float]]:
    module_genes = {gene for gene, score in affected_weights.items() if score > 0}
    direct_target_genes = set(target_weights)
    total_weight = sum(affected_weights.values()) or 1.0
    total_direct_weight = sum(target_weights.values()) or 1.0

    disease_frequency: dict[str, int] = defaultdict(int)
    for genes in disease_modules.values():
        for gene in genes:
            disease_frequency[gene] += 1
    disease_count = max(len(disease_modules), 1)
    gene_specificity = {
        gene: math.log((disease_count + 1.0) / (df + 0.5)) + 1.0
        for gene, df in disease_frequency.items()
    }
    weighted_denominator = sum(
        affected_weights.get(gene, 0.0) * gene_specificity.get(gene, 1.0)
        for gene in module_genes
    ) or total_weight

    raw_rows: list[dict[str, Any]] = []
    tdp_axis = axis_profile_from_genes(affected_weights)
    disease_axis_rows: list[dict[str, Any]] = []

    for disease, genes in disease_modules.items():
        overlap = sorted(module_genes.intersection(genes))
        direct_overlap = sorted(direct_target_genes.intersection(genes))
        specificity_overlap = (
            sum(affected_weights.get(gene, 0.0) * gene_specificity.get(gene, 1.0) for gene in overlap)
            / weighted_denominator
            * 100.0
        )
        direct_overlap_pct = (
            sum(target_weights.get(gene, 0.0) for gene in direct_overlap)
            / total_direct_weight
            * 100.0
        )
        p_value = hypergeom_tail(
            population_size=max(len(background), len(module_genes | genes)),
            success_population=len(genes),
            draw_count=len(module_genes),
            observed_success=len(overlap),
        )
        enrichment_score = min(-math.log10(max(p_value, 1e-300)) / 10.0 * 100.0, 100.0)
        disease_axis = axis_profile_from_genes(genes)
        tdp_match = cosine_similarity(tdp_axis, disease_axis)
        for axis, value in disease_axis.items():
            disease_axis_rows.append(
                {
                    "disease": disease,
                    "axis": axis,
                    "disease_axis_activity_0_100": round(value, 4),
                    "tdp_axis_activity_0_100": round(tdp_axis.get(axis, 0.0), 4),
                }
            )
        raw_rows.append(
            {
                "disease": disease,
                "disease_module_gene_count": len(genes),
                "overlap_count": len(overlap),
                "specificity_weighted_overlap_raw": specificity_overlap,
                "direct_target_overlap_raw": direct_overlap_pct,
                "enrichment_score": enrichment_score,
                "tdp_phenotype_match_score": tdp_match,
                "enrichment_p_value": p_value,
                "overlap_genes": ";".join(overlap[:80]),
                "direct_target_overlap_genes": ";".join(direct_overlap),
            }
        )

    max_specificity = max((row["specificity_weighted_overlap_raw"] for row in raw_rows), default=0.0)
    max_direct = max((row["direct_target_overlap_raw"] for row in raw_rows), default=0.0)
    for row in raw_rows:
        specificity_score = (
            row["specificity_weighted_overlap_raw"] / max_specificity * 100.0
            if max_specificity > 0
            else 0.0
        )
        direct_score = row["direct_target_overlap_raw"] / max_direct * 100.0 if max_direct > 0 else 0.0
        string_score = 0.70 * specificity_score + 0.15 * direct_score + 0.15 * row["enrichment_score"]
        final_score = 0.65 * string_score + 0.35 * row["tdp_phenotype_match_score"]
        row["string_module_disease_score"] = round(string_score, 4)
        row["final_indication_score"] = round(final_score, 4)
        row["specificity_weighted_overlap_score"] = round(specificity_score, 4)
        row["direct_target_overlap_score"] = round(direct_score, 4)
        row["specificity_weighted_overlap_raw"] = round(row["specificity_weighted_overlap_raw"], 6)
        row["direct_target_overlap_raw"] = round(row["direct_target_overlap_raw"], 6)
        row["tdp_phenotype_match_score"] = round(row["tdp_phenotype_match_score"], 4)
        row["enrichment_score"] = round(row["enrichment_score"], 4)
        row["enrichment_p_value"] = f"{row['enrichment_p_value']:.3e}"

    raw_rows.sort(key=lambda row: row["final_indication_score"], reverse=True)
    rows = raw_rows[:max_ranked_diseases]
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows, disease_axis_rows, tdp_axis


def write_markdown(path: Path, rows: list[dict[str, Any]]) -> None:
    lines = [
        "# Final PROTAC Indication Ranking",
        "",
        "This table was generated from a user-supplied target/off-target degradation profile, STRING network propagation, Enrichr disease modules, and automatic TDP phenotype-axis matching.",
        "",
        "| Rank | Disease module | Final | STRING module | TDP phenotype | Overlap genes |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['rank']} | {row['disease']} | {row['final_indication_score']:.2f} | "
            f"{row['string_module_disease_score']:.2f} | {row['tdp_phenotype_match_score']:.2f} | "
            f"{row['overlap_genes']} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-profile", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=root / "outputs" / "protac_indication_workflow")
    parser.add_argument("--limit-per-target", type=int, default=120)
    parser.add_argument("--max-candidate-disease-terms", type=int, default=800)
    parser.add_argument("--max-ranked-diseases", type=int, default=100)
    parser.add_argument(
        "--disease-keywords",
        nargs="+",
        default=DEFAULT_DISEASE_KEYWORDS,
        help="Disease-term keyword filter. Use 'all' for an unrestricted broad disease scan.",
    )
    parser.add_argument(
        "--enrichr-libraries",
        nargs="+",
        default=ENRICHR_LIBRARIES,
        help="Enrichr disease libraries to use. Defaults to the full configured library set.",
    )
    parser.add_argument("--include-oncology", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.out_dir / "cache"

    target_weights = load_targets(args.target_profile)
    target_payload = read_json(args.target_profile)
    exclude_oncology = bool(target_payload.get("exclude_oncology", True)) and not args.include_oncology

    string_edges = fetch_string_partners(
        target_weights,
        limit_per_target=args.limit_per_target,
        cache_path=cache_dir / "string_interaction_partners.json",
    )
    affected_weights, affected_rows = build_string_module(target_weights, string_edges, [])
    libraries = {library: fetch_enrichr_library(library, cache_dir) for library in args.enrichr_libraries}
    disease_modules, provenance_rows, background = build_candidate_modules(
        libraries,
        exclude_oncology=exclude_oncology,
        disease_keywords=args.disease_keywords,
        max_candidate_terms=args.max_candidate_disease_terms,
    )
    score_rows, disease_axis_rows, tdp_axis = score_modules(
        affected_weights=affected_weights,
        target_weights=target_weights,
        disease_modules=disease_modules,
        background=background,
        max_ranked_diseases=args.max_ranked_diseases,
    )

    write_csv(
        args.out_dir / "affected_protein_module.csv",
        affected_rows,
        ["protein", "string_module_weight", "module_role", "source_edges"],
    )
    write_csv(
        args.out_dir / "disease_module_terms.csv",
        provenance_rows,
        ["disease", "library", "matched_term", "term_gene_count", "new_genes_added"],
    )
    write_csv(
        args.out_dir / "tdp_axis_profile.csv",
        [
            {"axis": axis, "tdp_axis_activity_0_100": round(score, 4)}
            for axis, score in sorted(tdp_axis.items(), key=lambda item: item[1], reverse=True)
        ],
        ["axis", "tdp_axis_activity_0_100"],
    )
    write_csv(
        args.out_dir / "disease_axis_profiles.csv",
        disease_axis_rows,
        ["disease", "axis", "disease_axis_activity_0_100", "tdp_axis_activity_0_100"],
    )
    fields = [
        "rank",
        "disease",
        "final_indication_score",
        "string_module_disease_score",
        "tdp_phenotype_match_score",
        "specificity_weighted_overlap_score",
        "direct_target_overlap_score",
        "enrichment_score",
        "enrichment_p_value",
        "disease_module_gene_count",
        "overlap_count",
        "overlap_genes",
        "direct_target_overlap_genes",
    ]
    write_csv(args.out_dir / "final_indication_scores.csv", score_rows, fields)
    (args.out_dir / "final_indication_scores.json").write_text(
        json.dumps(score_rows, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_markdown(args.out_dir / "final_indication_scores.md", score_rows)

    print(
        json.dumps(
            {
                "ok": True,
                "targets": len(target_weights),
                "affected_module_proteins": len(affected_rows),
                "candidate_disease_modules": len(disease_modules),
                "ranked_diseases": len(score_rows),
                "top_disease": score_rows[0]["disease"] if score_rows else "",
                "out_dir": str(args.out_dir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
