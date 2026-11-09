#!/usr/bin/env python3
"""Download version-pinned discovery inputs; verify public archive and H5 hashes."""
import argparse
import hashlib
import json
import shutil
import urllib.request
import zipfile
from pathlib import Path

RECORD='https://zenodo.org/records/7760264/files/'
FILES={
    'DeconvolutionResults_ST_CRC_BelgianCohort.zip':'8a9604e51eff9f46a2d8fa4f56264980',
    'Pathology_SpotAnnotations.zip':'92882eb0db7c2cacc5077cb5fc9578fd'}
H5_SHA256='5c839f4e78b4f334e0c19dd3aef173c332b147b8418e854a0ffb6600e05da79a'

def digest(path, algorithm):
    h=hashlib.new(algorithm)
    with path.open('rb') as handle:
        while chunk:=handle.read(8*1024*1024):h.update(chunk)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    out=args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    archives=out/'archives'
    archives.mkdir(exist_ok=True)
    report=[]
    for name,expected_md5 in FILES.items():
        path=archives/name
        if not path.exists() or digest(path,'md5')!=expected_md5:
            temp=path.with_suffix('.partial')
            request=urllib.request.Request(RECORD+name+'?download=1',headers={'User-Agent':'COREAD-reproducibility/1.0'})
            with urllib.request.urlopen(request,timeout=90) as src,temp.open('wb') as dest:
                shutil.copyfileobj(src,dest,8*1024*1024)
            if digest(temp,'md5')!=expected_md5:
                raise ValueError(f'Public archive checksum mismatch: {name}')
            temp.replace(path)
        with zipfile.ZipFile(path) as z:
            for member in z.infolist():
                base=Path(member.filename).name
                if base in ['sp.h5ad','W_cell_density_q05.csv']:
                    destination=out/base
                elif base.startswith('Pathologist_Annotations_') and base.endswith('.csv'):
                    destination=out/'pathology'/base
                else:continue
                destination.parent.mkdir(exist_ok=True,parents=True)
                with z.open(member) as src,destination.open('wb') as dest:
                    shutil.copyfileobj(src,dest,8*1024*1024)
        report.append({'filename':name,'source_url':RECORD+name+'?download=1',
                       'md5':expected_md5,'sha256':digest(path,'sha256'),'bytes':path.stat().st_size})
        print(json.dumps(report[-1]),flush=True)
    actual=digest(out/'sp.h5ad','sha256')
    if actual!=H5_SHA256:raise ValueError('Extracted sp.h5ad differs from the analysis version')
    if len(list((out/'pathology').glob('*.csv')))!=14:raise ValueError('Expected 14 source pathology CSVs')
    (out/'source_download_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Verified source H5, 14 pathology tables and composition factors',flush=True)

if __name__=='__main__':main()
