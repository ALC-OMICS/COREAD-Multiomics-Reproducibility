#!/usr/bin/env python3
"""Unified entrypoint for dependency recomputation from archived processed data."""
from pathlib import Path
import argparse,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['processed','raw'],default='processed');p.add_argument('--output','--output-dir',dest='output',type=Path,required=True);a=p.parse_args()
    if a.mode=='raw':p.error('Original screen reads and fitted probability-training code were not supplied; use processed mode with the archived measurements.')
    out=a.output.resolve()
    if out.exists() and any(out.iterdir()):p.error('Choose a new or empty output directory')
    out.mkdir(parents=True,exist_ok=True);start=time.time();env=dict(os.environ,MPLBACKEND='Agg',MPLCONFIGDIR=str(out/'matplotlib_cache'))
    try:
        subprocess.run([sys.executable,str(ROOT/'run_all.py'),'--output',str(out/'numerical')],check=True,env=env)
        subprocess.run([sys.executable,str(ROOT/'scripts/export_current_panels.py'),'--numerical',str(out/'numerical'),'--output',str(out/'figures')],check=True,env=env)
        verification=json.loads((out/'numerical/verification_report.json').read_text());assert verification['status']=='PASS'
        status={'passed':True,'module':'01_dependency','mode':a.mode,'seconds':time.time()-start,'verification':verification,'limitations':['Original GO and cell-line release identifiers were not recovered.','Organoid probabilities use archived fitted coefficients; original fitting/training is not rerun.']}
    except Exception as exc:
        (out/'status.json').write_text(json.dumps({'passed':False,'module':'01_dependency','error':str(exc)},indent=2));raise
    (out/'status.json').write_text(json.dumps(status,indent=2));print(json.dumps(status,indent=2))
if __name__=='__main__':main()
