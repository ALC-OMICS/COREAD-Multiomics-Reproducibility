#!/usr/bin/env python3
"""Recompute single-cell inference from bundled pseudobulks or the public H5."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['processed','raw'],default='processed');p.add_argument('--output','--output-dir',dest='output',type=Path,required=True);p.add_argument('--h5',type=Path);p.add_argument('--expanded-dir',type=Path);a=p.parse_args()
    if a.mode=='raw' and not a.h5:p.error('Raw mode requires --h5 pointing to the public GSE178341 processed count file; see config/analysis_design.json for checksum.')
    out=a.output.resolve()
    if out.exists() and any(out.iterdir()):p.error('Choose a new or empty output directory')
    out.mkdir(parents=True,exist_ok=True);start=time.time();env=dict(os.environ,MPLBACKEND='Agg',MPLCONFIGDIR=str(out/'matplotlib_cache'))
    try:
        pseudobulk=ROOT/'input/pseudobulk'
        if a.mode=='raw':
            pseudobulk=out/'pseudobulk'
            subprocess.run([sys.executable,str(ROOT/'code/prepare_pseudobulk.py'),'--h5',str(a.h5.resolve()),'--output',str(pseudobulk),'--expanded-dir',str((a.expanded_dir or out/'expanded').resolve())],check=True,env=env)
        else:
            manifest=json.loads((pseudobulk/'manifest_sha256.json').read_text())
            for rel,digest in manifest.items():assert hashlib.sha256((pseudobulk/rel).read_bytes()).hexdigest()==digest,'Checkpoint checksum differs: '+rel
        subprocess.run([sys.executable,str(ROOT/'code/recompute_statistics.py'),'--pseudobulk',str(pseudobulk),'--output',str(out/'numerical')],check=True,env=env)
        subprocess.run([sys.executable,str(ROOT/'code/render_current_figure3.py'),'--results',str(out/'numerical'),'--output',str(out/'figures')],check=True,env=env)
        regression=json.loads((out/'numerical/numerical_regression.json').read_text());assert regression['passed']
        status={'passed':True,'module':'03_single_cell','mode':a.mode,'seconds':time.time()-start,'regression':regression,'limitations':['Inference code was reconstructed from source counts and the documented design; it is not the recovered original script.','Historical bootstrap random-stream scheduling is not established; newly computed intervals have separate labels.','The public processed matrix and source annotations are the starting point; alignment and source-study clustering are not reproduced.']}
    except Exception as exc:
        (out/'status.json').write_text(json.dumps({'passed':False,'module':'03_single_cell','error':str(exc)},indent=2));raise
    (out/'status.json').write_text(json.dumps(status,indent=2));print(json.dumps({'passed':True,'module':'03_single_cell','seconds':status['seconds'],'values_compared':regression['values_compared']},indent=2))
if __name__=='__main__':main()
