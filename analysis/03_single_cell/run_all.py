#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path

def run(cmd):
    print('+',' '.join(map(str,cmd)))
    subprocess.run(cmd,check=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--verify-only',action='store_true')
    ap.add_argument('--figures-only',action='store_true')
    ap.add_argument('--output',type=Path,default=Path('build'))
    args=ap.parse_args()
    root=Path(__file__).resolve().parent
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f'Output directory is not empty: {args.output}')
    args.output.mkdir(parents=True,exist_ok=True)
    if not args.figures_only:
        run([sys.executable,str(root/'code/verify_release.py'),'--root',str(root),'--report',str(args.output/'verification.txt')])
    if not args.verify_only:
        run([sys.executable,str(root/'code/make_figures.py'),'--root',str(root),'--output',str(args.output/'figures')])

if __name__=='__main__':main()
