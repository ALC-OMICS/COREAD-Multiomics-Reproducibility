#!/usr/bin/env python3
"""Read-only audit of the processed GSE178341 HDF5 and shipped metadata."""
from __future__ import annotations
import argparse, hashlib
from pathlib import Path
import h5py
import pandas as pd

EXPECTED_SHA = "f435bb2651ff5297d0c24a99daf58850ed67ae1ed6c5ef05fad48fa3f0186670"

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024), b""):
            h.update(b)
    return h.hexdigest()

def dec(a):
    return [x.decode() if isinstance(x,(bytes,bytearray)) else str(x) for x in a]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--h5", type=Path, required=True)
    ap.add_argument("--cluster", type=Path, required=True)
    ap.add_argument("--meta", type=Path, required=True)
    ap.add_argument("--mapping", type=Path, required=True)
    ap.add_argument("--skip-checksum", action="store_true")
    args=ap.parse_args()
    if not args.skip_checksum:
        d=sha256(args.h5)
        print("H5 SHA-256",d)
        assert d==EXPECTED_SHA, "H5 checksum mismatch"
    with h5py.File(args.h5,"r") as f:
        shape=tuple(int(x) for x in f["matrix/shape"][:])
        names=dec(f["matrix/features/name"][:])
        ids=dec(f["matrix/features/id"][:])
        nbc=len(f["matrix/barcodes"])
        print("H5 shape genes x cells:",shape)
        print("barcodes:",nbc)
    assert shape==(43113,370115)
    assert nbc==370115
    cluster=pd.read_csv(args.cluster)
    meta=pd.read_csv(args.meta)
    assert len(cluster)==370115 and len(meta)==370115
    assert cluster.sampleID.is_unique and meta.cellID.is_unique
    ann=cluster.merge(meta,left_on="sampleID",right_on="cellID",how="inner",validate="one_to_one")
    assert len(ann)==370115
    print("patients",ann.PID.nunique(),"tumor",ann.loc[ann.SPECIMEN_TYPE.eq('T'),'PID'].nunique(),"normal",ann.loc[ann.SPECIMEN_TYPE.eq('N'),'PID'].nunique())
    print("MMR tumor patients",ann.loc[ann.SPECIMEN_TYPE.eq('T'),['PID','MMRStatus']].drop_duplicates().MMRStatus.value_counts().to_dict())
    mp=pd.read_csv(args.mapping)
    assert len(mp)==466 and mp.Gene.nunique()==466
    for r in mp.itertuples(index=False):
        assert names[int(r.feature_row)]==r.H5_feature_name
        assert ids[int(r.feature_row)]==r.H5_feature_id
    print("feature mapping: PASS (466/466)")

if __name__=="__main__": main()
