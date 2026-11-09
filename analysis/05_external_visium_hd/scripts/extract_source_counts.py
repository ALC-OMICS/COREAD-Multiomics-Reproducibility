#!/usr/bin/env python3
"""Build compact reproducible checkpoints from verified public Visium HD H5 inputs."""
import argparse
import json
from pathlib import Path
import sys
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from counts import prepare_checkpoint

HASHES={
    'P1CRC':'e65dfe204f001364efa84b73c243c20605ab92aabb39e482afc1b0aadcfb6c02',
    'P2CRC':'445f02025d60a0aaa5bb2a0d9e4f5a476741d71a73c7d99c2397169cfba11b73',
    'P5CRC':'ee19eeeb56772ae6a83a06c8e7ca4d3ae00e3503de74680210fe6fe6e295023b'}

def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir',required=True,type=Path)
    parser.add_argument('--output-dir',type=Path,default=root/'data/counts')
    parser.add_argument('--sample',choices=list(HASHES),action='append')
    args=parser.parse_args()
    genes=pd.read_csv(root/'data/frozen_discovery_annotations.csv').Gene.tolist()
    for sample in args.sample or HASHES:
        qc=prepare_checkpoint(sample,args.raw_dir/f'{sample}_filtered_feature_bc_matrix.h5',
                              args.raw_dir/f'{sample}_Metadata.parquet.gz',genes,args.output_dir)
        if qc['source_h5_sha256'] != HASHES[sample]:
            raise ValueError(f'{sample}: source H5 differs from the verified analysis version')
        print(json.dumps(qc), flush=True)

if __name__=='__main__':main()
