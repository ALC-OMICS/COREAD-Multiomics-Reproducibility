# COREAD computational reproducibility

**Code/data-only release 0.2.1.** This repository contains executable analysis code,
processed measurement inputs, frozen numerical references and technical validation
records. It does **not** distribute the article manuscript, rendered figures,
figure legends, supplementary workbooks or journal submission documents.

The frozen candidate hierarchy is `Core47 -> Spatial43 -> Convergent9`, with
`Extended13` as a sensitivity set. The definitions are in
[`config/candidate_sets.json`](config/candidate_sets.json). No scientific selection
rule, random seed, statistical threshold or input measurement was changed by the
code-only packaging repair.

## Installation

Use Python 3.12 and the pinned dependencies. The earlier recovery environment was
Python 3.12.14; do not confuse that historical environment with the exact versions
recorded by a new validation run.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/verify_repository.py
python scripts/run_tests.py
```

On Windows, activate with `.venv\Scripts\activate`. Installation requires network
access; the default analyses use bundled processed inputs without downloading
source datasets. Optional source reconstruction has additional download and
resource requirements documented in each module.

## Execute and verify

Run from the repository root. Use a new output directory for each command.

```bash
# Small real-data demonstration (modules 04 and 07)
python run_all.py --demo --output-dir run_outputs/demo

# All seven numerical analysis modules
python run_all.py --output-dir run_outputs/full

# Refit spatial modules 04 and 05 from compact count inputs
python run_all.py --mode counts --analyses 04 05 --output-dir run_outputs/counts
```

Every run writes `run_report.json`. Successful execution requires a zero exit code,
`passed: true`, and successful module-specific checks. Generated graphics are local
run products in the chosen output directory; they are not versioned publication
figures and are excluded from this repository. Do not commit `run_outputs/`.

## Analytical modules

| Module | Computational scope |
|---|---|
| [`01_dependency`](analysis/01_dependency) | Cell-line and organoid measurements; recurrence and cross-system integration |
| [`02_genomics`](analysis/02_genomics) | Patient alteration endpoints, ranking, resampling and cohort comparisons |
| [`03_single_cell`](analysis/03_single_cell) | Patient-level pseudobulk inference and epithelial-composition sensitivity |
| [`04_discovery_spatial`](analysis/04_discovery_spatial) | Spatial statistics, patient aggregation and fixed candidate hierarchy |
| [`05_external_visium_hd`](analysis/05_external_visium_hd) | External spatial metrics, coverage and matched-background comparisons |
| [`06_crc16_niches`](analysis/06_crc16_niches) | Frozen host niches, bacterial associations and programme attribution |
| [`07_translational`](analysis/07_translational) | Checks of a frozen targeted evidence appraisal and deterministic mappings |

`source_data/` contains machine-readable numeric checkpoints and panel mappings,
not rendered figures or the article's supplementary workbooks. `results/` and
`reference/` inside modules contain fixed numerical regression targets. New
expected values are never generated from the rerun being assessed.

Five original cBioPortal OncoPrint XLSX exports remain in module 02's `input/`:
they are source measurement inputs, **not article supplementary tables**. Other
publication workbooks are neither included nor required. The former `ST11` CSV
checkpoint is retained byte-for-byte under the technical name
`discovery_spatial_expected.csv`; its function is numerical regression, not
publication formatting. The `modules/` directory retains legacy numeric
checkpoints; active entrypoints are in `analysis/`.

## Scope and provenance limits

Reproduction starts from published or archived **processed measurements**, not
from FASTQ alignment, de novo cell annotation or newly performed experiments.
Recovered and reconstructed implementations remain distinguished in each module.
The original GO query settings and cell-line release identifier were not
recoverable. Fitted organoid probability mappings are reconstructed, not retrained.
Some historical random-stream or unordered-set iteration details are unknown;
retained historical references are not presented as newly recovered algorithms.

Frozen numeric references are useful regression targets, not independent biological
validation. Cells and spatial bins do not increase patient replication. Missing
measurements remain missing; the external HMGB1 measurement is unavailable.
CRC16 microbial associations describe one specimen. The therapeutic appraisal is
frozen at its recorded dates and does not provide clinical treatment advice.

## Validation and maintenance

See [`validation/README.md`](validation/README.md) for the separation between
historical 0.2.0 records and current code-only checks. A report validates its stated
inputs, execution mode and environment, not all possible source-reconstruction
routes. Current reports should be read alongside the exact commit used.

`MANIFEST_SHA256.tsv` covers the current release, including the licence and
technical tests. `python scripts/verify_repository.py` checks hashes, required
inputs, candidate sets, numeric source mappings, Python syntax and the exclusion
of editorial attachments. It fails rather than silently accepting missing data.

After a reviewed maintenance change only, use
`python scripts/refresh_manifests.py` and repeat verification and tests. This
refreshes current inventories; historical manifest records and frozen numerical
reference values are not rewritten. Packaging uses `scripts/package_release.py`.

## Licence and citation

The existing [`LICENSE`](LICENSE) records MIT for original project code. Third-party
measurement exports and data retain their source terms; the software licence does
not replace them. See [`LICENSE_NOTICE.md`](LICENSE_NOTICE.md).

[`CITATION.cff`](CITATION.cff) identifies the software, maintainer and repository.
No article DOI or archival DOI is invented. Cite the exact tested release/commit
used for an analysis.
