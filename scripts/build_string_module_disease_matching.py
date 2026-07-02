#!/usr/bin/env python
"""STRING protein-module to disease-module matching.

This replaces OpenTargets or expression-matrix evidence in the degrader workflow with a
network-medicine style disease module matcher:

PROTAC target profile -> STRING affected protein module -> disease gene modules
from non-OpenTargets Enrichr libraries -> ranked autoimmune indications.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import requests


ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"
STRING_BASE = "https://string-db.org/api/json"
ENRICHR_LIBRARIES = [
    "Jensen_DISEASES_Curated_2025",
    "Jensen_DISEASES_Experimental_2025",
    "Jensen_DISEASES",
    "DisGeNET",
    "GWAS_Catalog_2025",
    "GWAS_Catalog_2023",
]

PATHWAY_AXIS_RULES = {
    "b cell receptor": {"b_cell_receptor_autoantibody": 1.0},
    "adaptive immune": {"b_cell_receptor_autoantibody": 0.4, "t_cell_receptor_signaling": 0.5, "myeloid_tlr_tnf_il6": 0.1},
    "fceri": {"mast_basophil_fceri": 1.0},
    "fc epsilon": {"mast_basophil_fceri": 1.0},
    "tslp": {"th2_cytokine": 1.0},
    "kit receptor": {"mast_basophil_fceri": 0.5, "epithelial_barrier_signaling": 0.2},
    "tyrosine kinase": {
        "b_cell_receptor_autoantibody": 0.3,
        "t_cell_receptor_signaling": 0.3,
        "mast_basophil_fceri": 0.2,
        "myeloid_tlr_tnf_il6": 0.2,
    },
    "intracellular signal transduction": {
        "b_cell_receptor_autoantibody": 0.25,
        "t_cell_receptor_signaling": 0.25,
        "mast_basophil_fceri": 0.2,
        "myeloid_tlr_tnf_il6": 0.2,
        "epithelial_barrier_signaling": 0.1,
    },
    "epidermal growth factor": {"epithelial_barrier_signaling": 0.8, "il23_il17_tissue_inflammation": 0.2},
    "epidermal": {"epithelial_barrier_signaling": 0.7, "il23_il17_tissue_inflammation": 0.3},
}


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def normalize_text(value: str) -> str:
    text = value.lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def load_targets(path: Path) -> dict[str, float]:
    payload = read_json(path)
    targets = {
        str(row.get("symbol", "")).upper(): safe_float(row.get("weight"))
        for row in payload.get("targets", [])
        if row.get("symbol") and not str(row.get("symbol", "")).upper().startswith("REPLACE_WITH")
    }
    if not targets:
        raise ValueError(
            "No usable target weights found. Provide a local target profile with real human gene symbols."
        )
    return targets


def fetch_string_partners(target_weights: dict[str, float], limit_per_target: int, cache_path: Path) -> list[dict[str, Any]]:
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))

    identifiers = "\r".join(target_weights)
    response = requests.post(
        f"{STRING_BASE}/interaction_partners",
        data={
            "identifiers": identifiers,
            "species": 9606,
            "limit": limit_per_target,
            "caller_identity": "codex-protac-workflow",
        },
        timeout=90,
    )
    response.raise_for_status()
    payload = response.json()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def build_string_module(
    target_weights: dict[str, float],
    string_edges: list[dict[str, Any]],
    local_network_rows: list[dict[str, str]],
) -> tuple[dict[str, float], list[dict[str, Any]]]:
    weights: dict[str, float] = {gene: weight * 100.0 for gene, weight in target_weights.items()}
    edge_notes: dict[str, list[str]] = defaultdict(list)

    for edge in string_edges:
        a = str(edge.get("preferredName_A", "")).upper()
        b = str(edge.get("preferredName_B", "")).upper()
        score = safe_float(edge.get("score"))
        for source, partner in [(a, b), (b, a)]:
            if source in target_weights and partner and partner != source:
                contribution = target_weights[source] * score * 100.0
                weights[partner] = max(weights.get(partner, 0.0), contribution)
                edge_notes[partner].append(f"{source}->{partner}:{target_weights[source]:.2f}*{score:.3f}")

    for row in local_network_rows:
        protein = row.get("protein", "").upper().strip()
        if not protein:
            continue
        score = safe_float(row.get("network_impact_score"))
        weights[protein] = max(weights.get(protein, 0.0), score)
        if row.get("source_edges"):
            edge_notes[protein].append(row["source_edges"])

    rows = []
    for protein, score in sorted(weights.items(), key=lambda item: item[1], reverse=True):
        rows.append(
            {
                "protein": protein,
                "string_module_weight": round(score, 4),
                "module_role": "direct_target" if protein in target_weights else "STRING_predicted_neighbor",
                "source_edges": " | ".join(edge_notes.get(protein, [])[:5]),
            }
        )
    return weights, rows


def fetch_enrichr_library(library: str, cache_dir: Path) -> dict[str, set[str]]:
    cache_path = cache_dir / f"{library}.txt"
    if not cache_path.exists():
        response = requests.get(
            f"{ENRICHR_BASE}/geneSetLibrary",
            params={"mode": "text", "libraryName": library},
            timeout=120,
        )
        response.raise_for_status()
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(response.text, encoding="utf-8")

    gene_sets: dict[str, set[str]] = {}
    for line in cache_path.read_text(encoding="utf-8").splitlines():
        parts = [part.strip() for part in line.split("\t")]
        if len(parts) < 3:
            continue
        term = parts[0]
        genes = {gene.upper() for gene in parts[2:] if re.fullmatch(r"[A-Za-z0-9_.-]+", gene)}
        if genes:
            gene_sets[term] = genes
    return gene_sets


def term_match_score(term: str, disease: str, disease_synonyms: dict[str, list[str]]) -> int:
    term_norm = normalize_text(term)
    synonyms = disease_synonyms.get(disease, [disease])
    best = 0
    for synonym in synonyms:
        syn_norm = normalize_text(synonym)
        if term_norm == syn_norm:
            best = max(best, 100)
        elif syn_norm in term_norm:
            best = max(best, 80 - min(25, len(term_norm.split()) - len(syn_norm.split())))
        elif term_norm in syn_norm and len(term_norm) >= 6:
            best = max(best, 55)
    return best


def build_disease_modules(
    diseases: list[str],
    disease_synonyms: dict[str, list[str]],
    libraries: dict[str, dict[str, set[str]]],
    max_terms_per_disease: int,
    max_genes_per_disease: int,
) -> tuple[dict[str, set[str]], list[dict[str, Any]], set[str]]:
    modules: dict[str, set[str]] = {}
    term_rows: list[dict[str, Any]] = []
    background: set[str] = set()
    for library_sets in libraries.values():
        for genes in library_sets.values():
            background.update(genes)

    for disease in diseases:
        candidates: list[tuple[int, int, str, str, set[str]]] = []
        for lib_index, (library, gene_sets) in enumerate(libraries.items()):
            for term, genes in gene_sets.items():
                score = term_match_score(term, disease, disease_synonyms)
                if score > 0:
                    candidates.append((score, -lib_index, library, term, genes))
        candidates.sort(key=lambda item: (item[0], item[1], -len(item[4])), reverse=True)

        selected = candidates[:max_terms_per_disease]
        module: set[str] = set()
        for score, _lib_priority, library, term, genes in selected:
            remaining_slots = max_genes_per_disease - len(module)
            if remaining_slots <= 0:
                break
            selected_genes = sorted(genes)[:remaining_slots]
            module.update(selected_genes)
            term_rows.append(
                {
                    "disease": disease,
                    "library": library,
                    "matched_term": term,
                    "term_match_score": score,
                    "term_gene_count": len(genes),
                    "genes_used_from_term": len(selected_genes),
                }
            )
        modules[disease] = module
    return modules, term_rows, background


def hypergeom_tail(population_size: int, success_population: int, draw_count: int, observed_success: int) -> float:
    if observed_success <= 0:
        return 1.0
    max_i = min(success_population, draw_count)
    log_den = math.lgamma(population_size + 1) - math.lgamma(draw_count + 1) - math.lgamma(population_size - draw_count + 1)
    probs = []
    for i in range(observed_success, max_i + 1):
        if draw_count - i > population_size - success_population:
            continue
        log_num = (
            math.lgamma(success_population + 1)
            - math.lgamma(i + 1)
            - math.lgamma(success_population - i + 1)
            + math.lgamma(population_size - success_population + 1)
            - math.lgamma(draw_count - i + 1)
            - math.lgamma(population_size - success_population - (draw_count - i) + 1)
        )
        probs.append(math.exp(log_num - log_den))
    return min(1.0, sum(probs))


def build_pathway_axis_scores(pathway_rows: list[dict[str, str]], disease_profiles: dict[str, dict[str, float]]) -> dict[str, float]:
    axis_scores: dict[str, float] = defaultdict(float)
    for row in pathway_rows:
        term = normalize_text(row.get("pathway_or_function", ""))
        priority = safe_float(row.get("priority_score"))
        for pattern, axis_weights in PATHWAY_AXIS_RULES.items():
            if pattern in term:
                for axis, weight in axis_weights.items():
                    axis_scores[axis] += priority * weight

    max_axis = max(axis_scores.values(), default=0.0)
    if max_axis > 0:
        axis_scores = {axis: value / max_axis * 100.0 for axis, value in axis_scores.items()}
    else:
        axis_scores = {}

    disease_scores: dict[str, float] = {}
    for disease, profile in disease_profiles.items():
        denom = sum(safe_float(v) for v in profile.values())
        if denom <= 0:
            disease_scores[disease] = 0.0
            continue
        disease_scores[disease] = sum(
            safe_float(need) * safe_float(axis_scores.get(axis))
            for axis, need in profile.items()
        ) / denom
    return disease_scores


def score_diseases(
    diseases: list[str],
    disease_modules: dict[str, set[str]],
    affected_weights: dict[str, float],
    target_weights: dict[str, float],
    background: set[str],
    pathway_scores: dict[str, float],
) -> list[dict[str, Any]]:
    module_genes = {gene for gene, score in affected_weights.items() if score > 0}
    total_weight = sum(affected_weights.values()) or 1.0
    direct_target_genes = set(target_weights)
    total_direct_weight = sum(target_weights.values()) or 1.0
    disease_count = max(len(diseases), 1)
    disease_frequency: dict[str, int] = defaultdict(int)
    for genes in disease_modules.values():
        for gene in genes:
            disease_frequency[gene] += 1
    gene_specificity = {
        gene: math.log((disease_count + 1.0) / (df + 0.5)) + 1.0
        for gene, df in disease_frequency.items()
    }
    weighted_module_denominator = sum(
        affected_weights.get(gene, 0.0) * gene_specificity.get(gene, 1.0)
        for gene in module_genes
    ) or total_weight
    raw_rows: list[dict[str, Any]] = []

    for disease in diseases:
        disease_genes = disease_modules.get(disease, set())
        overlap = sorted(module_genes.intersection(disease_genes))
        weighted_overlap_raw = sum(affected_weights.get(gene, 0.0) for gene in overlap) / total_weight * 100.0
        specificity_weighted_overlap_raw = (
            sum(
                affected_weights.get(gene, 0.0) * gene_specificity.get(gene, 1.0)
                for gene in overlap
            )
            / weighted_module_denominator
            * 100.0
        )
        direct_overlap_genes = sorted(direct_target_genes.intersection(disease_genes))
        direct_target_overlap_raw = (
            sum(target_weights.get(gene, 0.0) for gene in direct_overlap_genes)
            / total_direct_weight
            * 100.0
        )
        p_value = hypergeom_tail(
            population_size=max(len(background), len(module_genes | disease_genes)),
            success_population=len(disease_genes),
            draw_count=len(module_genes),
            observed_success=len(overlap),
        )
        enrichment_score = min(-math.log10(max(p_value, 1e-300)) / 10.0 * 100.0, 100.0)
        raw_rows.append(
            {
                "disease": disease,
                "disease_module_gene_count": len(disease_genes),
                "overlap_count": len(overlap),
                "weighted_overlap_raw_pct": weighted_overlap_raw,
                "specificity_weighted_overlap_raw_pct": specificity_weighted_overlap_raw,
                "direct_target_overlap_raw_pct": direct_target_overlap_raw,
                "enrichment_p_value": p_value,
                "enrichment_score": enrichment_score,
                "string_pathway_axis_score": pathway_scores.get(disease, 0.0),
                "overlap_genes": ";".join(overlap[:60]),
                "direct_target_overlap_genes": ";".join(direct_overlap_genes),
            }
        )

    max_overlap = max((row["weighted_overlap_raw_pct"] for row in raw_rows), default=0.0)
    max_specificity_overlap = max((row["specificity_weighted_overlap_raw_pct"] for row in raw_rows), default=0.0)
    max_direct_overlap = max((row["direct_target_overlap_raw_pct"] for row in raw_rows), default=0.0)
    for row in raw_rows:
        overlap_norm = row["weighted_overlap_raw_pct"] / max_overlap * 100.0 if max_overlap > 0 else 0.0
        specificity_overlap_norm = (
            row["specificity_weighted_overlap_raw_pct"] / max_specificity_overlap * 100.0
            if max_specificity_overlap > 0
            else 0.0
        )
        direct_overlap_norm = (
            row["direct_target_overlap_raw_pct"] / max_direct_overlap * 100.0
            if max_direct_overlap > 0
            else 0.0
        )
        final_score = (
            0.55 * specificity_overlap_norm
            + 0.15 * direct_overlap_norm
            + 0.10 * row["enrichment_score"]
            + 0.20 * row["string_pathway_axis_score"]
        )
        row["weighted_overlap_score"] = round(overlap_norm, 4)
        row["specificity_weighted_overlap_score"] = round(specificity_overlap_norm, 4)
        row["direct_target_overlap_score"] = round(direct_overlap_norm, 4)
        row["final_string_module_disease_score"] = round(final_score, 4)
        row["weighted_overlap_raw_pct"] = round(row["weighted_overlap_raw_pct"], 6)
        row["specificity_weighted_overlap_raw_pct"] = round(row["specificity_weighted_overlap_raw_pct"], 6)
        row["direct_target_overlap_raw_pct"] = round(row["direct_target_overlap_raw_pct"], 6)
        row["enrichment_p_value"] = f"{row['enrichment_p_value']:.3e}"
        row["enrichment_score"] = round(row["enrichment_score"], 4)
        row["string_pathway_axis_score"] = round(row["string_pathway_axis_score"], 4)

    raw_rows.sort(key=lambda row: row["final_string_module_disease_score"], reverse=True)
    for rank, row in enumerate(raw_rows, start=1):
        row["rank"] = rank
    return raw_rows


def write_markdown(path: Path, score_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# STRING Module Disease Matching",
        "",
        "This workflow uses STRING affected proteins/pathways plus non-OpenTargets disease gene modules.",
        "OpenTargets and expression-matrix disease evidence are not used.",
        "",
        "## Formula",
        "",
        "```text",
        "STRING affected module = direct PROTAC targets + STRING predicted neighbors, weighted by degradation depth and STRING confidence",
        "Disease module = matched disease gene sets from Jensen_DISEASES / DisGeNET / GWAS Catalog via Enrichr",
        "specificity weighted overlap = affected proteins overlapping disease module, downweighted when the same gene appears in many autoimmune disease modules",
        "direct target overlap = user-supplied degrader targets appearing directly in the disease module",
        "enrichment score = hypergeometric overlap significance",
        "STRING pathway axis score = STRING pathway impact projected to disease phenotype axes",
        "Final = 55% specificity weighted overlap + 15% direct target overlap + 10% enrichment + 20% STRING pathway axis",
        "```",
        "",
        "## Ranking",
        "",
        "| Rank | Disease | Final | Specific overlap | Direct target | Enrichment | STRING pathway axis | Overlap genes |",
        "|---:|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in score_rows:
        lines.append(
            f"| {row['rank']} | {row['disease']} | {row['final_string_module_disease_score']:.2f} | "
            f"{row['specificity_weighted_overlap_score']:.2f} | {row['direct_target_overlap_score']:.2f} | "
            f"{row['enrichment_score']:.2f} | "
            f"{row['string_pathway_axis_score']:.2f} | {row['overlap_genes']} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-profile", type=Path, default=root / "config" / "target_profile.template.json")
    parser.add_argument("--network-protein-impact", type=Path, default=root / "config" / "network_protein_impact.csv")
    parser.add_argument("--pathway-impact", type=Path, default=root / "config" / "pathway_impact.csv")
    parser.add_argument("--disease-profile", type=Path, default=root / "config" / "tdp_disease_profiles.template.json")
    parser.add_argument("--out-dir", type=Path, default=root / "outputs" / "string_module_disease_matching")
    parser.add_argument("--limit-per-target", type=int, default=120)
    parser.add_argument("--max-terms-per-disease", type=int, default=5)
    parser.add_argument("--max-genes-per-disease", type=int, default=500)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.out_dir / "cache"

    target_weights = load_targets(args.target_profile)
    local_network_rows = read_csv(args.network_protein_impact)
    pathway_rows = read_csv(args.pathway_impact)
    disease_profile_payload = read_json(args.disease_profile)
    disease_profiles = disease_profile_payload["disease_profiles"]
    disease_synonyms = disease_profile_payload.get("disease_synonyms", {})
    if not disease_profiles or any(str(name).startswith("REPLACE_WITH") for name in disease_profiles):
        raise ValueError(
            "No usable disease profiles found. Provide a local profile config with disease phenotype axes."
        )
    diseases = list(disease_profiles.keys())

    string_edges = fetch_string_partners(
        target_weights,
        limit_per_target=args.limit_per_target,
        cache_path=cache_dir / "string_interaction_partners.json",
    )
    affected_weights, affected_rows = build_string_module(target_weights, string_edges, local_network_rows)

    libraries = {
        library: fetch_enrichr_library(library, cache_dir)
        for library in ENRICHR_LIBRARIES
    }
    disease_modules, term_rows, background = build_disease_modules(
        diseases,
        disease_synonyms,
        libraries,
        max_terms_per_disease=args.max_terms_per_disease,
        max_genes_per_disease=args.max_genes_per_disease,
    )
    pathway_scores = build_pathway_axis_scores(pathway_rows, disease_profiles)
    score_rows = score_diseases(diseases, disease_modules, affected_weights, target_weights, background, pathway_scores)

    write_csv(
        args.out_dir / "string_affected_protein_module.csv",
        affected_rows,
        ["protein", "string_module_weight", "module_role", "source_edges"],
    )
    write_csv(
        args.out_dir / "disease_module_terms.csv",
        term_rows,
        ["disease", "library", "matched_term", "term_match_score", "term_gene_count", "genes_used_from_term"],
    )
    write_csv(
        args.out_dir / "string_module_disease_scores.csv",
        score_rows,
        [
            "rank",
            "disease",
            "final_string_module_disease_score",
            "specificity_weighted_overlap_score",
            "specificity_weighted_overlap_raw_pct",
            "direct_target_overlap_score",
            "direct_target_overlap_raw_pct",
            "direct_target_overlap_genes",
            "weighted_overlap_score",
            "weighted_overlap_raw_pct",
            "enrichment_score",
            "enrichment_p_value",
            "string_pathway_axis_score",
            "disease_module_gene_count",
            "overlap_count",
            "overlap_genes",
        ],
    )
    (args.out_dir / "string_module_disease_scores.json").write_text(
        json.dumps(score_rows, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_markdown(args.out_dir / "string_module_disease_matching.md", score_rows)

    print(
        json.dumps(
            {
                "ok": True,
                "diseases": len(score_rows),
                "affected_module_proteins": len(affected_rows),
                "disease_module_terms": len(term_rows),
                "top_disease": score_rows[0]["disease"] if score_rows else "",
                "out_dir": str(args.out_dir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
