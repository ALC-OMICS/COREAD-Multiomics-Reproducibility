> Historical 0.2.0 documentation, not current execution instructions.
> Publication attachments mentioned below are not distributed in 0.2.1.

# COREAD dependency atlas

**Cell lines · organoids · cross-system evidence integration**  
Reproducibility package assembled on **3 October 2026** for Andrés López-Cortés.

This repository contains the current **1,720-symbol analysis**, using the corrected
**1,165-model cell-line matrix**, the recovered individual NextGen organoid Chronos
estimates, final supplementary workbooks, source data, quantitative SVG figures,
code, model metadata, and computational audit records. The obsolete aggregate-bound
analysis and earlier graphics are isolated under `archive/`.

## Results at a glance

| Item | Current result |
|---|---:|
| Input gene symbols | 1,720 |
| Supplied cell-line models | 1,165 |
| COREAD cell-line models | 63 |
| Primary other-cancer comparators | 1,099 |
| Source-labelled Non-Cancerous models, excluded from primary comparator | 3 |
| Genes evaluable in cell lines | 1,620 |
| Cell-line candidates | 335 |
| Higher COREAD call frequency, Fisher BH q < 0.05 | 23 |
| Primary COREAD organoid models, 3DO | 16 |
| Expanded COREAD organoid models, 3DO + 2DO | 18 |
| Genes evaluable in organoids | 1,636 |
| Primary organoid candidates | 324 |
| Same-model Chronos / re-estimated-probability concordance | 188 |
| Genes evaluable in both systems | 1,613 |
| Candidates in both systems | 193 |
| Candidates in the full observed union | 466 |

These are **exploratory candidate rules**, not FDR-controlled dependency discoveries.
Comparative FDR describes between-group differences, not the error rate of the
335-, 324-, or 466-gene catalogues.

## Quick start

No account credentials, web downloads, commercial tools or spreadsheet APIs are
needed to replay the included data. Python **3.13** is recommended; the included
package versions were tested with Python 3.13.5.

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell alternative: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python run_all.py --figures
```

Outputs go to `build/`; frozen reference workbooks and source data are not replaced.
Use `python run_all.py --output /path/to/a/separate/folder --figures` for another
output location. Runtime and memory depend on the computer: the full analysis
includes 10,000 lineage-bootstrap draws and 2,000 integration-bootstrap draws.

To check the distributed snapshots without rerunning the numerical analyses:

```bash
python run_all.py --verify-only
python -m unittest discover -s tests -v
```

The entrypoint runs the organoid reconstruction, cell-line statistics and integration,
optional figures, and comparison to the final spreadsheet cell values. Reproduction
uses numerical tolerances; XLSX styling and ZIP timestamps are not regenerated.

## What can be reproduced

| Component | Included implementation and scope |
|---|---|
| Cell-line statistics | Exact recovered script, repository-relative paths; Fisher/BH, CIs, conditional odds ratios, lineage comparisons, bootstrap and related-pair sensitivities |
| Organoid effects and candidate selection | Recomputed from archived individual gene-effect values; checked against both recovered query and genome-wide extractions |
| Organoid probabilities | Numerically reconstructed from archived per-screen logistic coefficients; **not a new KDE/EM fit** |
| Organoid comparative context | Welch/BH, Mann–Whitney/BH and format-adjusted OLS-HC3/BH recomputed from recovered 108-screen data |
| Integration | Candidate union/intersection, missingness-aware evidence patterns, recurrence sensitivities, Jaccard/Dice, rank association and model-bootstrap intervals |
| Figures | Four editable, vector-only SVGs; source data and plotting scripts included |
| Supplementary Excel files | Frozen byte-identical final files; their cached gene-table values are verified against recalculated CSV outputs |

**The original empirical-Bayesian KDE/EM training and five-fold internal-validation
script was not recovered.** Its recorded parameters, 648 training controls,
18,435 × 18 genome-wide inputs, fitted coefficients and reported internal metrics
are retained. The replay reconstructs probabilities from the fitted coefficients
and does not claim to rerun those original fits or their held-out checks. No raw
CRISPR read-count processing or Chronos refitting is performed.

## Data and analysis rules

Both model systems are assessed in parallel from the original query list. A
cell-line candidate requires a supplied dependency probability **>0.5 in ≥2 COREAD
models**. An organoid candidate requires an unmodified Chronos score **<−0.5 in ≥2
of the 16 3DO models**. These are different outcomes and different recurrence
floors; their probabilities, effects and P values are never pooled.

Cell-line denominators exclude missing measurements per gene. The main comparison
excludes only the three explicit `Non-Cancerous` source labels; the sensitivity
includes all 1,102 non-COREAD models. Source diagnoses are preserved, not silently
re-curated. Per-lineage BH correction covers 38,755 eligible gene-by-lineage tests.

In the integration, 462 candidates are jointly evaluable: 193 shared + 139
cell-line-only + 130 organoid-only. The full union is **466**, because another three
cell-line candidates lack organoid data and one organoid candidate lacks cell-line
data. Missing data are not classified as a negative dependency result.

## Repository map

```text
coread-dependency-atlas/
├── README.md / LEEME.md       Guidance in English and Spanish
├── run_all.py                Reproduction entrypoint
├── requirements.txt          Tested package versions
├── CITATION.cff              Repository citation metadata
├── LICENSING.md              Source / code licensing boundaries
├── data/
│   ├── inputs/               Supplied source matrix and audit inputs
│   ├── organoids/            Recovered scores and re-estimated probabilities
│   ├── metadata/             Query genes, model cohorts, EB coefficients/controls
│   └── reference/            Recovered upstream metadata and older context output
├── scripts/                  Numerical reconstruction and plotting code
├── results/
│   ├── tables/               Three final, single-sheet Excel workbooks
│   ├── csv/                  Machine-readable results
│   └── gene_sets/            Plain-text candidate lists
├── figures/                  Canonical SVGs, PNG previews, source data
├── docs/                     Methods, results, provenance notes and references
├── provenance/               SHA-256 manifest and verification reports
├── tests/                    Automated consistency tests
└── archive/                  Superseded analyses and original code/figures
```

## Important limitations

The supplied cell-line release and paired Chronos matrix remain **unidentified /
unavailable**. No version was inferred from the number of models. Organoid
probabilities are re-estimates, not the unavailable published probability matrix.
Model-ID overlap is checked, but donor independence and complete model relatedness
are not established. The two documented related COREAD pairs are only a limited
sensitivity analysis. None of the analyses demonstrates microbiome causality,
clinical efficacy, normal-tissue selectivity or independent biological validation.

Axis membership was supplied upstream. Gene-by-axis annotation and BIC/microbiome,
mutation, drug or treatment datasets are **not part of this dependency repository**.

## Sources and publication

See recorded references (`docs/REFERENCES.md`, historical path), reproducibility scope (`docs/REPRODUCIBILITY.md`, historical path),
methods (`docs/METHODS.md`, historical path), results (`docs/RESULTS.md`, historical path), and licensing notice (`LICENSING.md`, historical path).
No upstream data or code is relabelled as exclusively owned by this project.
The cell-line input lacks a release and redistribution licence declaration;
resolve those metadata before representing the data as an openly licensed release.

For GitHub, extract the distribution ZIP and upload the **contents** of this folder,
so that this README is at the repository root. Do not upload only the outer ZIP.
No repository has been created or published automatically.
