#!/usr/bin/env python3
"""Reproduce numerical results, optionally regenerate SVGs, then verify frozen outputs.

No web access, authentication, external uploads, or source-overwriting operations
are performed. Results are written under --output (default: build/).
"""
from pathlib import Path
import argparse,os,subprocess,sys,json
ROOT=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'build')
    p.add_argument('--figures',action='store_true',help='Also generate the four SVG/PNG figures')
    p.add_argument('--verify-only',action='store_true',help='Check included snapshots without rerunning numerical statistics')
    args=p.parse_args();out=args.output.resolve()
    if args.verify_only:
        subprocess.run([sys.executable,str(ROOT/'scripts/verify_results.py')],check=True,cwd=ROOT);return
    if out==ROOT or out in [ROOT/'data',ROOT/'results']:p.error('Choose a separate output directory')
    out.mkdir(parents=True,exist_ok=True);logs=out/'logs';logs.mkdir(exist_ok=True)
    env=dict(os.environ,COREAD_OUTPUT_DIR=str(out),MPLBACKEND='Agg',COREAD_FIGURE_DIR=str(out/'figures'),OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
    stages=['reproduce_organoids.py','reproduce_cell_lines.py']
    if args.figures:stages+=['plot_organoids.py','plot_cell_lines_and_integration.py']
    for stage in stages:
        print(f'Running {stage}',flush=True)
        with (logs/(stage+'.log')).open('w') as f:
            completed=subprocess.run([sys.executable,str(ROOT/'scripts'/stage)],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        if completed.returncode:
            print((logs/(stage+'.log')).read_text()[-6000:]);raise SystemExit(completed.returncode)
    subprocess.run([sys.executable,str(ROOT/'scripts/verify_results.py'),'--build',str(out)],cwd=ROOT,env=env,check=True)
    print('Complete. Numerical outputs and logs: '+str(out))
if __name__=='__main__':main()
