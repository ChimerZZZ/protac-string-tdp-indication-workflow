# STRING-TDP Network Medicine Workflow for PROTAC Indication Prioritization

This repository implements a reproducible computational workflow for prioritizing autoimmune and inflammatory indications from a user-supplied targeted protein degradation profile.

The workflow is designed for early-stage translational research, where the key question is:

> Given a quantitative degrader perturbation profile, which disease modules are most consistent with the predicted proteomic and phenotypic effects of the molecule?

The public README intentionally does **not** disclose project-specific input values, target-level measurements, disease scores, or ranking results. Those should be kept in local configuration and output files.

## Scope

Targeted protein degradation differs from conventional inhibition because a degrader can produce graded and multi-target effects. A useful workflow should therefore account for:

- target degradation depth across intended targets and off-targets,
- network-level propagation from direct targets to affected proteins,
- disease module proximity or overlap,
- pathway-level immune context,
- phenotype-level match between the degrader effect and disease biology.

This repository focuses on the workflow itself. It does not ship project-specific input data, result tables, manuscript figures, or proprietary biological conclusions.

## Model Stack

The final workflow combines three public/open-source components:

- [STRING API](https://string-db.org/help/api/) for protein interaction partners and pathway enrichment.
- [Enrichr](https://maayanlab.cloud/Enrichr/) for non-OpenTargets disease gene-set libraries, including Jensen Diseases, DisGeNET, and GWAS Catalog resources.
- [NetMedPy](https://github.com/menicgiulia/NetMedPy)-style network medicine logic for module proximity and overlap scoring.

The workflow does not require expression matrices, OpenTargets disease evidence, therapeutic-precedence evidence, or a compound-specific training set. Users provide their own degrader target profile and disease phenotype profile.

## Workflow Summary

```text
User-supplied target degradation profile
        |
        v
STRING affected protein module
        |
        v
Disease modules from Enrichr libraries
        |
        v
Network medicine module matching
        |
        v
Target degradation phenotype specificity gate
        |
        v
Final indication prioritization table
```

## Input Interface

The workflow expects a local JSON configuration describing a target degradation profile. This file should contain target identifiers and normalized degradation weights.

Example schema:

```json
{
  "targets": [
    {
      "symbol": "TARGET_A",
      "weight": 1.0,
      "note": "local evidence note"
    }
  ],
  "exclude_oncology": true
}
```

Project-specific target identities, assay values, and disease scores should remain in local configuration and output files.

Use the files in `config/*.template.json` as schemas only. For a real run, create local files such as:

```text
config/target_profile.local.json
config/tdp_disease_profiles.local.json
```

Local profile files should not be committed to a public repository.

## Model Details

### 1. STRING Affected Protein Module

Direct degraded targets are expanded through STRING interaction partners. The affected protein module is weighted by target degradation strength and STRING interaction confidence.

```text
STRING module weight(gene)
= max over source targets [target degradation weight x STRING interaction score x 100]
```

### 2. Disease Modules

Disease modules are assembled from Enrichr disease gene-set libraries. The current workflow supports Jensen Diseases, DisGeNET, and GWAS Catalog-derived disease modules.

The mapping step is configured through synonym and disease-term matching rules inside the scripts and local configuration files.
Candidate diseases and disease synonyms should be supplied in the local TDP/disease profile file.

### 3. Specificity-Weighted Disease Module Matching

The affected STRING module is matched against each disease module using:

- specificity-weighted overlap,
- direct target overlap,
- hypergeometric enrichment,
- STRING pathway axis score.

Generic immune genes that appear across many disease modules are downweighted using an IDF-like disease-frequency penalty:

```text
gene specificity(g)
= log((N_diseases + 1) / (disease_frequency(g) + 0.5)) + 1
```

The STRING disease-module score is:

```text
STRING disease-module score
= 55% specificity-weighted overlap
 + 15% direct target overlap
 + 10% enrichment score
 + 20% STRING pathway axis score
```

### 4. Target Degradation Phenotype Matching

The target degradation profile is also mapped onto disease-relevant phenotype axes such as immunoreceptor signaling, lymphocyte biology, innate immune signaling, cytokine response, tissue inflammation, and epithelial/barrier biology.

```text
TDP axis activity
= sum(target degradation weight x target-axis effect prior)
```

```text
TDP phenotype match
= weighted coverage of disease phenotype axes by TDP axis activity
```

### 5. Final Score

The final score combines protein-module matching with phenotype specificity:

```text
Final score
= 65% STRING disease-module score
 + 35% TDP phenotype specificity score
```

This design keeps STRING-based proteome-to-disease matching as the primary model while using TDP phenotype matching as a specificity gate.

## Reproducible Commands

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the one-command workflow with a local target/off-target profile:

```bash
python scripts/run_protac_indication_workflow.py \
  --target-profile config/target_profile.local.json
```

This minimal mode automatically builds disease modules from Enrichr libraries and infers TDP phenotype-axis matching from the STRING-expanded affected protein module.
By default, candidate disease terms are filtered toward autoimmune and inflammatory indications and oncology terms are excluded. Use `--disease-keywords all` and `--include-oncology` for a broader disease scan.

Run the advanced workflow with an explicit local disease/TDP profile:

```bash
python scripts/build_tdp_phenotype_matching.py \
  --target-profile config/target_profile.local.json \
  --profile-config config/tdp_disease_profiles.local.json

python scripts/build_string_module_disease_matching.py \
  --target-profile config/target_profile.local.json \
  --disease-profile config/tdp_disease_profiles.local.json

python scripts/build_string_tdp_final_scores.py
```

The first STRING/Enrichr run may take several minutes because public gene-set libraries and STRING interaction results are cached locally.

## Output Interface

The workflow writes local outputs under `outputs/`. Typical outputs include:

- affected protein module tables,
- disease module provenance tables,
- STRING disease-module matching scores,
- TDP phenotype matching scores,
- final indication prioritization tables.

Generated outputs are intentionally excluded from version control by default.

## Repository Structure

```text
config/
  target_profile.template.json
  tdp_disease_profiles.template.json

scripts/
  build_tdp_phenotype_matching.py
  build_string_module_disease_matching.py
  build_string_tdp_final_scores.py
  run_protac_indication_workflow.py
```

## Limitations

This is a research-prioritization workflow, not a clinical efficacy predictor.

- STRING interactions are predicted or curated network evidence, not measured degrader-induced proteomics.
- Disease modules can contain generic immune genes; specificity weighting reduces but does not eliminate this issue.
- TDP phenotype priors should be updated when measured phosphoproteomics, transcriptomics, cytokine, or cell-state data become available.
- Final scores are relative prioritization scores and should be interpreted as hypothesis-generating outputs.

## Citation / Attribution

If adapting this workflow, cite the upstream resources used for protein network inference, disease gene-set construction, and network medicine analysis:

- STRING database / STRING API
- Enrichr and the underlying disease gene-set libraries
- NetMedPy and the network medicine proximity framework
