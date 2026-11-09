#!/usr/bin/env python3
"""Run all module and repository unit tests in isolated subprocesses."""
from pathlib import Path
import argparse, json, os, subprocess, sys, time
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report',type=Path)
    a=p.parse_args()
    folders=sorted((ROOT/'analysis').glob('*/tests'))
    if (ROOT/'tests').is_dir(): folders.append(ROOT/'tests')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MPLBACKEND='Agg')
    results=[]
    for folder in folders:
        started=time.perf_counter()
        result=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(folder),'-p','test_*.py','-v'],cwd=ROOT,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        entry={'suite':folder.relative_to(ROOT).as_posix(),'passed':result.returncode==0,'returncode':result.returncode,'seconds':round(time.perf_counter()-started,3),'output':result.stdout}
        results.append(entry)
        print(('PASS ' if entry['passed'] else 'FAIL ')+entry['suite'],flush=True)
        if not entry['passed']: print(result.stdout,flush=True)
    report={'passed':bool(results) and all(r['passed'] for r in results),'suites':results,'python':sys.version}
    if a.report:
        a.report.parent.mkdir(parents=True,exist_ok=True)
        a.report.write_text(json.dumps(report,indent=2)+'\n')
    return 0 if report['passed'] else 1
if __name__=='__main__':sys.exit(main())
