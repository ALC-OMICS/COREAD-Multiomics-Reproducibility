# Reproducibility

## Source object
- File: `ST.seu.rds`
- SHA-256: `bb04b5f057f173f789fd9a36c2b30b17f31329b677e8812dfe19014aa43babc3`
- Tissue spots analysed: 3,106

The source RDS is not redistributed in this GitHub-ready folder.

## Primary design
- Host niche definition: host SCT expression only.
- PCA components used for graph construction: 20.
- k-nearest neighbours: 20.
- Primary Louvain resolution: 0.8.
- Primary host niches: 13.
- Spatial35: 35 frozen genes.
- Convergent6: MYC, FASN, SOX9, BCL2L1, IRS2, CTNNB1.
- Primary bacterial family: 20 taxa × 13 niches = 260 tests.
- Multiple-testing correction: Benjamini–Hochberg within the prespecified family.
- GSE178341: source annotations retained; 62 tumour patients; patient-level pseudobulk contrasts.

## Software environment used for downstream analysis
- Python 3.13.5
- NumPy 2.3.5
- pandas 2.2.3
- SciPy 1.17.0
- scikit-learn 1.8.0
- statsmodels 0.14.6
- NetworkX 3.6.1
- h5py 3.15.1

## Important interpretation boundary
CRC_16 is a within-tumour spatial discovery analysis. The spatial results do not constitute patient-level external validation and are not direct evidence of clinical checkpoint-inhibitor response.
