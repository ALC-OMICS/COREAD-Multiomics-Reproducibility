# Module 01 — cell-line and organoid dependency analysis

Code/data-only interface: no supplementary workbook is required or shipped.
The cell-line probabilities and organoid individual-effect matrices remain the
calculation inputs. Previously archived numerical CSV files in `results/csv/`
are fixed regression references, not output newly generated as its own expected
answer. Stored field vocabulary and organoid annotations are reused where stated;
cell-line statistics and organoid candidate/contextual statistics are recomputed.

From the repository root:

```bash
python analysis/01_dependency/run.py --mode processed --output-dir run_outputs/dependency
python analysis/01_dependency/scripts/verify_results.py
python -m unittest discover -s analysis/01_dependency/tests -v
```

`run.py` executes numerical reconstruction, checks fixed candidate counts, compares
866,880 record values across the three primary references, and produces `status.json`.
Read `verification_report.json` for scope. Additional organoid contextual values
and reconstructed fitted probability mappings are checked separately. Plotting
helpers write only generated local outputs under the selected run directory.

The fixed universe/order is in `data/metadata/query_genes_1720.txt`. Model identities,
metadata, fitted coefficients and controls are separate files under `data/metadata/`.
`data/inputs/Texto pegado(1).txt` preserves the full archived probability matrix;
its release identifier and original GO query settings were not recovered.
The older 250-model matrix is retained solely for recorded overlap checks.

Organoid probabilities reconstruct the archived fitted mapping; the unavailable
KDE/EM training procedure is not refitted or represented as recovered original
code. The two model systems use different measurements. Candidate recurrence is
not itself evidence of CRC preference or normal-tissue selectivity.

The retained `archive/` and `provenance/` records identify historical code and
source hashes. Their old filenames can refer to publication products no longer
included. Active checksums are in `provenance/sha256_manifest.json` and the root
manifest. The default pipeline does not execute historical workbook helpers.
