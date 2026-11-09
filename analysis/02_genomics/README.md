# Module 02 — patient alteration analysis

This module reproduces patient-level endpoints and their numerical regression
checks without requiring the article's supplementary workbooks or a `tables/`
folder. The five XLSX files under `input/` are original cBioPortal OncoPrint
measurement exports, not manuscript supplements; their cached source values are
read without reevaluating formulas.

From the repository root:

```bash
python analysis/02_genomics/run.py --mode processed --output-dir run_outputs/genomics
python -m unittest discover -s analysis/02_genomics/tests -v
```

The wrapper creates a separate working directory, parses source exports, reruns
the patient analyses, verifies permutation probabilities and BH families, and
compares six CSV files plus nine NPZ files against the fixed numerical references.
A new or empty output directory is required. Generated graphics are local run
products, never required packaged publication figures.

`input/union_candidates_466.txt` fixes the candidate universe. Patient coverage,
mutation definitions and modality-specific denominators follow the archived input
annotations; they are not silently harmonised across cohorts. Missing driver
labels remain unavailable. Paired PanCancer comparisons are not independent
patient replication. The expanded AC-ICAM acquisition route is optional and
records its pinned public-data provenance; the default route uses bundled
processed summaries and the frozen patient-call extract.

`results/` contains immutable numerical regression targets. Tests check cohort
coverage, unique-patient unions, missingness, empirical probabilities, rank
bootstrap arrays and numeric display-source consistency. They do not require
an Excel export of the article tables. Historical manifests and reports under
`provenance/` are distinguished from the current root checksum inventory.
