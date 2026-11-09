# Raw HDF5 input

The 1.203-GB processed expression matrix is intentionally **not committed** because it exceeds GitHub's normal per-file limit.

Expected file:

`GSE178341_crc10x_full_c295v4_submit.h5`

Expected SHA-256:

`f435bb2651ff5297d0c24a99daf58850ed67ae1ed6c5ef05fad48fa3f0186670`

Obtain the processed matrix from GEO accession **GSE178341** (Pelka et al., Cell 2021), then place it in this directory. Alternatively, if the file is stored as sequential binary parts named `*.part001`, `*.part002`, ... use `python code/prepare_h5.py --parts-dir data/raw`.

The processed cluster and patient metadata are included in `input/` as gzip-compressed CSV files because they are small enough for GitHub and are required to reproduce the published annotations used in this analysis.
