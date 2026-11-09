# Changes in the GitHub-ready packaging

## 1.0.0 — 2026-10-06

- Recovered the original 94-entry patient-analysis archive and preserved its original tables and numerical results.
- Fixed the input path in `parse_exports.py` to use the repository input folder rather than an authoring-environment absolute path. Original code bytes retained in `archive/original_analysis_code/`.
- Added an isolated-output runner, numerical regression comparison and read-only tests.
- Added reconstructed code and data for the later compact composition plots. Source statistics and the primary 30-gene order are unchanged; headings, margins and endpoint explanations are clarified.
- Added a source-exact compact permutation chart with all displayed P/FDR values, including weaker results, and an explicit PanCancer robustness label. Retained the same 15 displayed contrasts and original FDR families.
- Added documentation of endpoint differences, normalization, assay availability, bootstrap/permutation scope and optional authoring-specific Excel export.
- Excluded generative-AI quantitative chart recreations, unused font files and execution caches.
- No new genes were selected, no alteration/driver definitions harmonized retrospectively, and no new primary analysis performed by changing scientific parameters.
