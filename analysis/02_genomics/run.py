#!/usr/bin/env python3
"""Unified genomics entrypoint: original patient-export inference and figures."""
from pathlib import Path
import argparse,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['processed','raw'],default='processed');ap.add_argument('--output','--output-dir',dest='output',type=Path,required=True);a=ap.parse_args()
    if a.mode=='raw':ap.error('This module starts from archived patient exports and pinned external summaries. Optional AC-ICAM acquisition is documented in code/fetch_acicam.py; original variant calling is outside this repository.')
    out=a.output.resolve()
    if out.exists() and any(out.iterdir()):ap.error('Choose a new or empty output directory')
    out.mkdir(parents=True,exist_ok=True);start=time.time();env=dict(os.environ,MPLBACKEND='Agg',MPLCONFIGDIR=str(out/'matplotlib_cache'))
    try:
        subprocess.run([sys.executable,str(ROOT/'run_all.py'),'--output',str(out/'numerical')],check=True,env=env)
        subprocess.run([sys.executable,str(ROOT/'code/export_current_panels.py'),'--numerical',str(out/'numerical'),'--output',str(out/'figures')],check=True,env=env)
        regression=json.loads((out/'numerical/provenance/numerical_regression.json').read_text());summary=json.loads((out/'numerical/results/analysis_summary.json').read_text());assert regression['result']=='PASS'
        status={'passed':True,'module':'02_genomics','mode':a.mode,'seconds':time.time()-start,'regression':regression,'summary':summary,'limitations':['Original patient exports and archived AC-ICAM summaries are processed inputs.','PanCancer is a paired TCGA representation check, not an independent validation cohort.','AC-ICAM broad coding/splice mutation calls do not provide driver-specific annotations.']}
    except Exception as exc:
        (out/'status.json').write_text(json.dumps({'passed':False,'module':'02_genomics','error':str(exc)},indent=2));raise
    (out/'status.json').write_text(json.dumps(status,indent=2));print(json.dumps({'passed':True,'module':'02_genomics','seconds':status['seconds'],'numerical_files_checked':len(regression['checks'])},indent=2))
if __name__=='__main__':main()
