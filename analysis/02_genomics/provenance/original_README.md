> Historical 0.2.0 documentation, not current execution instructions.
> Publication attachments mentioned below are not distributed in 0.2.1.

# COREAD patient alteration atlas

Reproducible patient-level molecular profiling of **466 frozen dependency candidates** in colorectal adenocarcinoma, with external cohort comparisons and permutation-based concordance assessment.

This repository covers the **patient-analysis stage** (Supplementary Tables 7–9), not a new DepMap screen or a rerun of the preceding cell-line/organoid analysis. The upstream functional evidence and the 466-gene input set are retained for traceability. This is an exploratory research resource, not a clinical test or a list of proven therapeutic targets.

**Guía en español (`LEEME.md`, historical path)** · **Methods (`docs/METHODS.md`, historical path)** · **Results (`docs/RESULTS.md`, historical path)** · **Figure guide (`docs/FIGURES.md`, historical path)** · **Reproducibility (`docs/REPRODUCIBILITY.md`, historical path)** · **Sources and citations (`CITATION.md`, historical path)**

## Design

| Role | Studies | Available patient-analysis material |
|---|---|---|
| Primary | TCGA COAD-GDC + READ-GDC | 631 unique patients; 572 jointly mutation/CNA-profiled (421 colon, 151 rectum). |
| External colon cohort | AC-ICAM (`coad_silu_2022`) | 281 DNA patients, 456 exact candidate symbols; 348 RNA patients, 451 exact symbols. Archived source summaries, raw-file hashes and top-30 individual call extract are supplied. |
| External colon cohort | CPTAC-2 prospective (`coad_cptac_2019`) | 110 patients in the export; 105 jointly mutation/CNA-profiled. |
| Analytical robustness only | TCGA PanCancer Atlas | Shared TCGA patients; paired coverage and tests. **Not independent validation.** |

The principal DNA rank is the proportion of unique jointly profiled patients with a **putative-driver mutation OR amplification OR deep deletion**. Frequencies are not added. AC-ICAM has no driver annotations, so its broad coding/splice mutation endpoint is explicitly separate. RNA remains descriptive because the source thresholds/references were not harmonized; protein frequencies were not calculated without verified gene-level coverage.

The three statistical families contain **18 primary external tests**, **18 colon-only sensitivity tests** and **10 paired-TCGA robustness tests**. Of 50 planned contrasts, 46 are evaluable and each has **10,000 stored null draws**. Four unavailable AC-ICAM driver-specific contrasts have no P value. Ranking stability uses **2,000 patient bootstrap draws**. These quantities assess concordance or stability, not gene-level clinical validation.

## Start here

Read-only verification of the shipped data:

```bash
python -m pip install -r requirements.txt
python run_all.py --verify-only
```

Complete offline numerical reanalysis and figure generation:

```bash
python run_all.py
```

The runner creates `build/` and leaves the reference `input/`, `results/`, `tables/` and `figures/` unchanged. If that directory already contains files, use a different destination:

```bash
python run_all.py --output build-second-run
```

To regenerate figures only, without recomputing inferential statistics:

```bash
python run_all.py --figures-only --output build-figures
```

The numerical path needs Python, NumPy, SciPy and Matplotlib. Exact library versions used for the bundled reproduction check are pinned in `requirements.txt`; Python 3.13.5 was used. No notebook, API key or network request is required once dependencies are installed. The default pipeline does **not** require the authoring-specific Excel runtime.

### What is reproduced offline?

The five categorical Excel source exports are reparsed. From those calls and the **archived AC-ICAM source-derived summaries**, the code recalculates frequencies, Wilson intervals, Fisher/McNemar comparisons, BH corrections, all 46 permutation tests and the 2,000 stratified patient-bootstrap draws. The top-30 AC-ICAM patient call extract is provided separately for its OncoPrint.

This is **not** raw-sequencing reprocessing, a new driver annotation or a fresh download of the AC-ICAM genome-wide files. The optional `code/fetch_acicam.py` retrieves pinned official files and checks their hashes; this network-dependent route is separate and was not rerun when this repository was packaged.

### Excel files

The three one-sheet workbooks are supplied **unchanged**. The runner verifies their values against the reproduced CSV outputs; it does not need to rewrite them. `code/build_tables.py` preserves the original authoring/export procedure, which needs `artifact_tool` and is optional. Portable numerical reproduction does not depend on this package. There is no fabricated replacement of missing values.

## Files

| Directory | Contents |
|---|---|
| `input/` | Five original patient-level exports; frozen 466 genes; upstream functional evidence CSV and original input symbols. |
| `downloads/` | Recovered AC-ICAM count summaries, RNA/cytoband annotations, top-30 individual calls and transfer-verification artifacts. Not the full raw AC-ICAM sequencing/matrix collection. |
| `code/` | Original numerical analysis with a documented input-path portability fix; source acquisition, tables, verification and reconstructed compact-plot code. |
| `results/` | Complete primary/external tables, patient-call arrays, all endpoint summaries, all null draws, bootstrap ranks/frequencies and the display-only top-30 list. |
| `tables/` | Supplementary Tables **7, 8 and 9**, each with one sheet. |
| `figures/` | Seven original numerical figures plus new compact profiles and a compact permutation summary, with PNG previews and CSV source data. |
| `docs/` | Methods, results, source/endpoint dictionary, figure semantics, limitations and reproduction guide. |
| `provenance/` | Input audits, source hashes, environment, original checks and current repository validation reports. |
| `tests/` | Read-only numerical/data-integrity tests. |
| `archive/` | Byte-preserved original scripts/provenance and earlier delivered numerical SVG versions. **Not a source of alternative statistical results.** |

## Reading the figures correctly

The conventional OncoPrints show patients × genes. The newer stacked bars show a **unit-sum normalization of five marginal alteration frequencies**; they are not conventional OncoPrints and a full-width bar does not mean that every patient is altered. The percentages to the right are separate, unnormalized DNA-union patient frequencies. Mini-tiles indicate the presence of an alteration category, **not the five biological axes**.

For TCGA and CPTAC-2 the right percentage uses the primary putative-driver/CNA union. For AC-ICAM it uses the **broader coding/splice/CNA union**. No protein or structural-variant segments are invented. RNA categories are not tumor–normal differential-expression results.

The compact permutation summary shows 15 selected gene-label tests, retaining the P and BH values from their **original complete families**. Its broad DNA row is not the driver/CNA endpoint of the principal ranking. All other tests remain in Table 9 and the full numerical results. Null intervals represent the central permutation range, **not confidence intervals for the observed correlation**.

Only numerical, data-linked charts are included in `figures/`. Generative-AI recreations of quantitative charts are excluded because they do not preserve the underlying measurements.

## Interpretation and release conditions

All **466 genes remain in the candidate catalogue**. The top 30 are a display choice, not a newly validated subset. A significant Fisher result indicates cohort heterogeneity; a nonsignificant result does not prove equivalence. Global permutation significance is not an external-validation P value for an individual gene. See limitations (`docs/LIMITATIONS.md`, historical path).

The source cohorts and archived public study identifiers are preserved. No public GitHub repository has been created by this packaging step. Before public release, the research team should review attribution, consent/data-use restrictions applicable to the source resources and the desired license. This package does not assign a blanket open license to third-party datasets. See rights notice (`LICENSE_NOTICE.md`, historical path).
