# Validation scope — code/data-only release 0.2.1

`archive/0.2.0/` contains unchanged records from the original recovery package.
Their timings, dependency installation and deliverable checks do not certify
this code-only version. Paths in those records are historical, not current links.

Current validation reports, when present in `current/`, record the tests actually
executed after the repair: clean environment details, all seven processed-input
modules, spatial count-based refits and the unit-test suites. A PASS requires
successful module checks, not merely a zero process exit code.

Run these checks yourself from the root:

```bash
python scripts/verify_repository.py
python scripts/run_tests.py --report run_outputs/unit_tests.json
python run_all.py --output-dir run_outputs/processed
python run_all.py --mode counts --analyses 04 05 --output-dir run_outputs/counts
```

The dependency module compares fresh calculations with existing frozen CSV
references. The genomic module compares six CSV files and nine NPZ files and
independently verifies permutation probabilities/BH adjustments. No publication
workbook is required by either analysis. The same source measurement arrays,
statistical thresholds, seeds and candidate sets are retained.

The count-based route refits spatial statistics from bundled compact counts;
it does not repeat raw-read alignment or source-study annotation. The full
GSE178341 H5 and CRC16 RDS reconstruction paths require optional external source
files and are not covered merely by a processed/count-mode PASS.

Generated plots remain local run outputs. No rendered graphics, supplementary
workbooks or manuscript text are uploaded with validation reports. Historical
provenance files keep their recorded filenames/checksums even when those old
publication attachments are no longer distributed.
