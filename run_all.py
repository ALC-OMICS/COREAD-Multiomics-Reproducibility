#!/usr/bin/env python3
"""Run the seven analysis modules against bundled processed inputs."""
from __future__ import annotations
import argparse, importlib.metadata, json, os, platform, subprocess, sys, tempfile, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MODULES = {
    '01': '01_dependency',
    '02': '02_genomics',
    '03': '03_single_cell',
    '04': '04_discovery_spatial',
    '05': '05_external_visium_hd',
    '06': '06_crc16_niches',
    '07': '07_translational',
}

def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix='.'+path.name+'.')
    with os.fdopen(fd, 'w') as f:
        json.dump(value, f, indent=2); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(tmp,path)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--analyses', nargs='+', choices=list(MODULES), help='Module numbers; default: all seven')
    ap.add_argument('--demo', action='store_true', help='Run modules 04 and 07 as a small real-data demonstration')
    ap.add_argument('--mode', choices=['processed','counts'], default='processed', help='counts refits modules 04/05 from bundled compact counts; other modules retain their processed-input contracts')
    ap.add_argument('--output-dir', type=Path, default=ROOT/'run_outputs')
    args = ap.parse_args()
    chosen = args.analyses or (['04','07'] if args.demo else list(MODULES))
    out = args.output_dir.resolve()
    if out == ROOT or (out.is_relative_to(ROOT) and not out.is_relative_to(ROOT/'run_outputs')):
        ap.error('Use run_outputs/ or a directory outside the repository; source folders are read-only.')
    out.mkdir(parents=True, exist_ok=True); (out/'logs').mkdir(exist_ok=True)
    env=os.environ.copy()
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
        env[key]='1'
    env['MPLBACKEND']='Agg'; env['PYTHONHASHSEED']='0'
    env.setdefault('MPLCONFIGDIR',str(out/'matplotlib_cache'))
    versions={}
    for name in ['numpy','scipy','pandas','matplotlib','h5py','statsmodels','scikit-learn','pyarrow','openpyxl','patsy']:
        try: versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: versions[name]=None
    report={'release':'0.2.1','demo':args.demo,'python':sys.version,'platform':platform.platform(),
            'dependencies':versions,'mode':args.mode,'modules':[],'passed':False}
    start=time.perf_counter()
    for key in chosen:
        module=MODULES[key]; target=out/module
        if target.exists() and any(target.iterdir()):
            raise SystemExit(f'Use a new output directory; {target} is non-empty.')
        module_mode='counts' if args.mode=='counts' and key in ['04','05'] else 'processed'
        command=[sys.executable,str(ROOT/'analysis'/module/'run.py'),'--mode',module_mode,'--output-dir',str(target)]
        print(f'Running {module}...',flush=True); t=time.perf_counter()
        with (out/'logs'/(module+'.log')).open('w') as log:
            completed=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        status_path=target/'status.json'
        status=json.loads(status_path.read_text()) if status_path.exists() else {'passed':False,'error':'Module did not produce status.json'}
        entry={'module':module,'seconds':round(time.perf_counter()-t,3),'returncode':completed.returncode,
               'passed':completed.returncode==0 and status.get('passed') is True,
               'status':status,'log':'logs/'+module+'.log'}
        report['modules'].append(entry); report['seconds']=round(time.perf_counter()-start,3)
        report['passed']=all(x['passed'] for x in report['modules']) and len(report['modules'])==len(chosen)
        atomic_json(out/'run_report.json',report)
        print(f"  {'PASS' if entry['passed'] else 'FAIL'} ({entry['seconds']} s)",flush=True)
        if not entry['passed']:
            print(f'See {out/entry["log"]}',file=sys.stderr);return 1
    print(f'Completed {len(chosen)} modules in {report["seconds"]} s. Report: {out/"run_report.json"}')
    return 0

if __name__=='__main__': sys.exit(main())
