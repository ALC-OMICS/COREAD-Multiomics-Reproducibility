> Historical 0.2.0 documentation, not current execution instructions.
> Publication attachments mentioned below are not distributed in 0.2.1.

# COREAD single-cell dependency prioritization

Patient-level single-cell transcriptomic localization and prioritization of **466 colorectal adenocarcinoma (COREAD) dependency candidates** using the processed **GSE178341** atlas from Pelka et al. (Cell, 2021).

This repository covers the single-cell stage that follows the cell-line/organoid dependency screen and the patient genomic alteration analysis. It contains the numerical results underlying **Supplementary Table 10 and Figure 3A–D**, together with source metadata, gene matching, analysis design, figure source data and reproducibility checks.

**Guía en español (`LEEME.md`, historical path)** · **Methods (`docs/METHODS.md`, historical path)** · **Results (`docs/RESULTS.md`, historical path)** · **Figures (`docs/FIGURES.md`, historical path)** · **Reproducibility (`docs/REPRODUCIBILITY.md`, historical path)** · **References (`docs/REFERENCES.md`, historical path)**

## Main design

| Component | Design |
|---|---|
| Candidate universe | 466 genes fixed before single-cell analysis |
| Source atlas | GSE178341, Pelka et al., Cell 2021 |
| Processed cells | 370,115 annotated cells |
| Tumor patients | 62: 34 MMRd and 28 MMRp |
| Adjacent-normal comparison | 36 patient-matched epithelial pairs |
| Primary inference | Patient-level pseudobulk, not cell-level pseudoreplication |
| Primary test | Two-sided Wilcoxon signed-rank, BH across 466 genes |
| Cell-count sensitivity | >=100 epithelial cells/specimen, 35 pairs |
| Composition sensitivity | >=20 cells/shared epithelial state, >=2 states/patient, 33 pairs |
| Tumor epithelium vs microenvironment | 62 matched patients |
| MMR comparison | MMRd n=34 vs MMRp n=28, two-sided Mann-Whitney U |

The source annotations comprise seven top-level lineages: epithelial, T/NK/ILC, myeloid, B, plasma, stromal and mast cells. Tumor-derived epithelial cells are **not automatically equivalent to independently confirmed malignant cells**; no copy-number-based malignant-cell reclassification was imposed.

## Principal results

All 466 candidates were represented in the expression matrix: **455 exact symbol matches and 11 stable-Ensembl recoveries**. In the paired tumor-derived versus adjacent-normal epithelial analysis, 334 genes had BH FDR <0.05: 219 with positive median effects, 112 with negative effects and 3 with zero median despite significant rank-distribution differences. The >=100-cell sensitivity identified 220 positive signals, and the epithelial-composition-balanced sensitivity identified 198 positive signals.

The matched tumor epithelium versus aggregate non-epithelial tumor comparison identified 206 positive and 179 negative FDR-significant genes. Dominant-lineage assignments across the complete catalogue were stromal (140), epithelial (113), T/NK/ILC (76), myeloid (65), mast (41), B (19) and plasma (12).

Integration of four prespecified evidence dimensions generated **65 Tier A**, **134 Tier B** and **267 Tier C** genes. Requiring Tier A membership plus a positive composition-balanced tumor-versus-normal epithelial effect at FDR <0.05 defined the **47-gene spatial-follow-up core**. The 47-gene set is a stringent evidence-integration subset for subsequent spatial analysis, **not a list of validated therapeutic targets**.

## Start here

Install the numerical environment:

```bash
python -m pip install -r requirements.txt
```

Verify the shipped release:

```bash
python run_all.py --verify-only
```

Regenerate Figure 3A–D from the released source CSV files:

```bash
python run_all.py --figures-only --output build-figures
```

Run both verification and figure regeneration:

```bash
python run_all.py --output build
```

The reference `results/`, `tables/` and `figures/` directories are never overwritten by the runner.

## Raw HDF5 source

The processed GSE178341 HDF5 matrix is 1.203 GB and is **not committed** because it exceeds GitHub's standard per-file limit. Its expected SHA-256 is:

```text
f435bb2651ff5297d0c24a99daf58850ed67ae1ed6c5ef05fad48fa3f0186670
```

Place the source matrix at `data/raw/GSE178341_crc10x_full_c295v4_submit.h5` for local source inspection. If you have sequential binary parts, `code/prepare_h5.py` reconstructs the file and verifies the digest. The two processed annotation CSVs used in the analysis are supplied in `input/` in gzip-compressed form.

The repository is a reproducibility package for the **reported processed-data analysis**; it does not re-run FASTQ alignment, the original source-study quality-control/integration pipeline, or the original clustering and cell annotation.

## Repository structure

| Directory | Contents |
|---|---|
| `input/` | Frozen 466-gene list, upstream functional/patient evidence, exact feature matching, compressed GSE178341 annotations |
| `results/` | Complete single-cell gene-level results, sensitivities, lineage and epithelial-state summaries, Tier A and Core 47 lists |
| `tables/` | Plain one-sheet Supplementary Table 10 |
| `figures/` | Figure 3A–D in SVG/PNG plus numerical source CSVs |
| `config/` | Frozen analysis design and thresholds |
| `code/` | HDF5 reconstruction/source inspection, release validation and figure regeneration |
| `docs/` | Methods, Results, figure semantics, data dictionary, limitations, references and reproducibility notes |
| `provenance/` | Checksums, environment, raw-source manifest and validation reports |
| `tests/` | Read-only integrity tests |

## Interpretation

The evidence layers are intentionally kept separate. Transcript abundance does not demonstrate functional dependency in an individual patient, and dependency does not imply a therapeutic window. MMR-associated expression is secondary and is not part of Tier A/Core 47 membership. The single-cell priority rules were exploratory and not preregistered. Multiple-testing corrections apply within the explicitly defined 466-gene comparison families, not to a unified study-wide FDR.

No public GitHub repository is created by this package. Before publication, review the desired code license and the source-data terms applicable to third-party resources; see LICENSE_NOTICE.md (`LICENSE_NOTICE.md`, historical path).
