# Citation and source attribution

This folder is an unpublished research-analysis snapshot (repository package version 1.0.0, 2026-10-06). No public repository URL, archived release DOI or complete authorship list has been assigned here. Add the final repository URL, release DOI and confirmed author list before citing the software as a published resource. Cite the underlying cohorts and statistical resources separately.

## Cohort and portal sources used by the analysis

- cBioPortal COAD-GDC: https://www.cbioportal.org/study/summary?id=coad_tcga_gdc
- cBioPortal READ-GDC: https://www.cbioportal.org/study/summary?id=read_tcga_gdc
- AC-ICAM pinned DataHub source: https://github.com/cBioPortal/datahub/tree/86690e1ed9752b1dcd50b5657f5f05eafa4b6b78/public/coad_silu_2022
- CPTAC-2: https://www.cbioportal.org/study/summary?id=coad_cptac_2019
- PanCancer representation: https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
- HGNC cytoband source recorded in the archive: https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

Actual byte-level inputs and raw-source hashes are identified in `provenance/`, not inferred from a current website version. The source URLs do not establish an unrecorded original export date or query threshold.

## References from the analysis documentation

1. Cerami E, Gao J, Dogrusoz U, et al. The cBio Cancer Genomics Portal: an open platform for exploring multidimensional cancer genomics data. Cancer Discovery. 2012;2:401–404. https://doi.org/10.1158/2159-8290.CD-12-0095
2. Gao J, Aksoy BA, Dogrusoz U, et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Science Signaling. 2013;6:pl1. https://doi.org/10.1126/scisignal.2004088
3. The Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. Nature. 2012;487:330–337. https://doi.org/10.1038/nature11252
4. Roelands J, Kuppen PJK, Ahmed EI, et al. An integrated tumor, immune and microbiome atlas of colon cancer. Nature Medicine. 2023;29:1273–1286. https://doi.org/10.1038/s41591-023-02324-5
5. Vasaikar S, Huang C, Wang X, et al. Proteogenomic analysis of human colon cancer reveals new therapeutic opportunities. Cell. 2019;177:1035–1049.e19. https://doi.org/10.1016/j.cell.2019.03.030
6. Hoadley KA, Yau C, Hinoue T, et al. Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. Cell. 2018;173:291–304.e6. https://doi.org/10.1016/j.cell.2018.03.022
7. Mermel CH, Schumacher SE, Hill B, et al. GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biology. 2011;12:R41. https://doi.org/10.1186/gb-2011-12-4-r41
8. Wilson EB. Probable inference, the law of succession, and statistical inference. Journal of the American Statistical Association. 1927;22:209–212. https://doi.org/10.1080/01621459.1927.10502953
9. Newcombe RG. Interval estimation for the difference between independent proportions: comparison of eleven methods. Statistics in Medicine. 1998;17:873–890. https://doi.org/10.1002/(SICI)1097-0258(19980430)17:8%3C873::AID-SIM779%3E3.0.CO;2-I
10. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society, Series B. 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
11. McNemar Q. Note on the sampling error of the difference between correlated proportions or percentages. Psychometrika. 1947;12:153–157. https://doi.org/10.1007/BF02295996
12. Phipson B, Smyth GK. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. Statistical Applications in Genetics and Molecular Biology. 2010;9:Article 39. https://doi.org/10.2202/1544-6115.1585
13. Seal RL, Braschi B, Gray K, et al. Genenames.org: the HGNC resources in 2023. Nucleic Acids Research. 2023;51:D1003–D1009. https://doi.org/10.1093/nar/gkac888
14. Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods. 2020;17:261–272. https://doi.org/10.1038/s41592-019-0686-2

## Required interpretation

TCGA PanCancer is not an independent validation cohort. AC-ICAM mutation calls in this analysis are not driver-annotated. The analysis does not calculate harmonized RNA or protein validation statistics. Global permutation concordance is not individual-gene validation. The source-specific uncertainties and original processing versions are documented in the archive.
