# Provenance records

Files inherited from the original patient-analysis archive (including
`final_verification.json`, `environment.json`, `input_audit.json` and
`sha256_manifest.json`) document that original authoring snapshot. Their original
copies are also identifiable in `archive/original_analysis_provenance/`.

The original `sha256_manifest.json` is **not** the checksum manifest for the new
GitHub package after documented path and plotting-code changes. Use the root
`MANIFEST_SHA256.tsv` for the current package.

New records prefixed `repository_` describe the packaging run, test environment,
clean numerical reanalysis, regression and current tests. The clean-run report
explicitly identifies AC-ICAM as supplied summaries/top-30 calls; it does not
claim a fresh download or raw-sequencing rerun. An additional margin-only change
to the compact permutation plot was rerendered after the numeric reanalysis;
this does not change the test values.

`acicam_public_source.json` records original public-file LFS hashes and the pinned
source revision. Local recovered CSV summaries are different files with their
own hashes in the package manifest.
