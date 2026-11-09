# Source Data for Figures 1–7

`panel_index.csv` maps all 31 panels to their source tables, descriptions, units and missingness rules. Each figure subdirectory contains ordinary CSV tables suitable for independent analysis. Figure 1 has panels A–E; Figures 2–5 and 7 have A–D; Figure 6 has A–F.

`SourceData_Figure_1.tsv` through `SourceData_Figure_7.tsv` are lossless containers for convenient figure-level inspection. Each container has explicit record types: FORMAT, FIGURE, PANEL, TABLE, HEADER, DATA, END_TABLE and END_PANEL. A TABLE record gives the original relative filename, its SHA-256 and row count. HEADER and DATA records preserve every original string cell, including blanks and NA markers. They are not transformed numerical summaries. Use the individual CSV files when a single rectangular data table is needed.

`source_data_manifest.json` records the source tables, row counts, file hashes and successful cell-by-cell round-trip checks. Rebuild separately with:

```bash
python scripts/build_source_data.py --output-dir run_outputs/source_data
```

Candidate membership and availability are retained in the underlying data. For Figure 4A, all 466 source genes are included even though only 462 have finite plotting coordinates. For Figure 5, the nine Convergent9 genes are listed in each tumour and the three HMGB1 observations remain missing. Heatmap masks, set denominators and categorical zero values have the meanings specified in the panel index and legends.

The source tables support numerical reproduction of the figures. Original and regenerated layouts are distinguished in `figures`, and individual module outputs contain additional statistics beyond the displayed panels.
