# Executed validation — release 0.2.0

Validation date: 31 October 2026. The central run used a fresh Python 3.12.14 virtual environment without system-site-packages.

## Central execution

- Clean dependency installation: **23.351 seconds**, successful.
- Complete processed-input run: **7/7 modules passed**, **77.91 seconds**.
- Real-data demo: **2/2 modules passed**, **6.214 seconds**.

These are observed times on the recorded environment. They are not a hardware-independent performance claim. See the JSON reports and per-module logs for the executed checks.

| Module | Observed time (seconds) | Status |
|---|---:|---|
| 01_dependency | 34.598 | PASS |
| 02_genomics | 20.922 | PASS |
| 03_single_cell | 5.593 | PASS |
| 04_discovery_spatial | 4.957 | PASS |
| 05_external_visium_hd | 5.809 | PASS |
| 06_crc16_niches | 4.945 | PASS |
| 07_translational | 1.072 | PASS |

## Numerical evidence

| Analysis | Independent checks |
|---|---|
| Dependency | 866,880 recalculated values compared with three archived tables; 335 cell-line, 324 organoid, 193 shared and 466 union candidates. |
| Genomics | 15 numerical reference files, including 46 × 10,000 permutation draws and the 2,000-resample bootstrap arrays. |
| Single-cell | 64 fields × 466 genes = 29,824 comparisons; maximum difference approximately 4.55e-13; Core47 recovered exactly. |
| Discovery spatial | All 159 non-key ST11 fields × 466 genes reproduced. Source H5 and compact-count refits additionally verify normalisation, composition adjustment and spatial statistics. Historical null vectors remain distinct from new canonical-order draws. |
| External spatial | All current ST12 sections checked; full count/library normalisation, 1,595,565 bins and 40,000 freshly seeded matched-null draws. |
| CRC16 | 559 Spatial43 gene/niche effects, 91 immune programme effects, 260 planned primary tests, 240 evaluable revised CMH tests and six exact C9 patient-lineage contrasts. |
| Translational | Ten records with classes 3A/6B/1D; 32 coverage cells comprising 8 direct, 5 conditional and 19 zero values. |

## Deliverable integrity

All 14 supplementary worksheets were exported and checked cell by cell after a CSV round trip (24,618 rows). All 31 figure panels map to present source tables. Seven figure-level TSV containers preserve 51 table blocks with exact string-cell round trips. Python sources compile under the tested environment. ZIP CRC checks and member-level SHA-256 comparisons are performed by the packaging script; extraction from the smaller parts is tested separately.

## Interpretation

A PASS records the stated comparison; it does not imply recovery of an unavailable historical script, clinical validation of a treatment or author approval of publication declarations. See `docs/PROVENANCE_AND_LIMITS.md` and the module records for the explicit boundaries.

## Additional count-based central run

Modules 04 and 05 passed the central `--mode counts` route in 77.828 seconds, using bundled counts and full-library metadata without network access. See `counts_run_report.json`. Module 04 refits the composition model; module 05 recomputes all gene-level spatial metrics before set-level inference.
