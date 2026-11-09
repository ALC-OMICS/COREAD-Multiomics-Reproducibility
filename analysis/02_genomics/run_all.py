#!/usr/bin/env python3
"""Offline reproduction in an isolated output directory; reference files are immutable.

Examples:
    python run_all.py --verify-only
    python run_all.py
    python run_all.py --output /absolute/path/new-output
    python run_all.py --figures-only --output build-figures
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def run_step(work: Path, name: str) -> dict:
    start = time.monotonic()
    log = work / 'provenance' / (Path(name).stem + '.log')
    env = {**os.environ, 'MPLBACKEND': 'Agg', 'OPENBLAS_NUM_THREADS': '1',
           'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}
    with log.open('w', encoding='utf-8') as handle:
        cp = subprocess.run([sys.executable, str(work / 'code' / name)],
                            cwd=work, env=env, stdout=handle, stderr=subprocess.STDOUT)
    if cp.returncode:
        raise RuntimeError(f'{name} failed (exit {cp.returncode}).\n{log.read_text()[-8000:]}')
    result = {'script': name, 'returncode': 0, 'seconds': round(time.monotonic()-start, 3)}
    print(f"OK  {name} ({result['seconds']} s)", flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'build',
                        help='New or empty output directory; never overwrites reference outputs.')
    parser.add_argument('--verify-only', action='store_true', help='Run archived-data unit checks; no recomputation.')
    parser.add_argument('--figures-only', action='store_true', help='Rebuild plots from supplied numerical results.')
    args = parser.parse_args()
    if args.verify_only and args.figures_only:
        parser.error('Choose only one of --verify-only and --figures-only.')
    if args.verify_only:
        subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT/'tests'), '-v'],
                       cwd=ROOT, check=True)
        return
    work = args.output.expanduser().resolve()
    if work == ROOT or ROOT.is_relative_to(work):
        parser.error('Output cannot be the repository root or one of its parents.')
    if work.exists() and any(work.iterdir()):
        parser.error(f'{work} is not empty. Choose a new directory; no existing files were removed.')
    work.mkdir(parents=True, exist_ok=True)
    for sub in ['code', 'input', 'downloads', 'config']:
        shutil.copytree(ROOT/sub, work/sub, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__', 'raw_acicam'))
    for sub in ['results', 'figures', 'provenance']:
        (work/sub).mkdir(exist_ok=True)
    shutil.copy2(ROOT/'provenance/acicam_public_source.json', work/'provenance/acicam_public_source.json')
    if args.figures_only:
        shutil.copytree(ROOT/'results', work/'results', dirs_exist_ok=True)
        stages = []
    else:
        # AC-ICAM all-gene inputs are archived counts; the 30-gene individual call
        # extract is also provided for the corresponding patient OncoPrint.
        shutil.copy2(ROOT/'results/AC_ICAM_top30_calls.npz', work/'results/AC_ICAM_top30_calls.npz')
        stages = ['parse_exports.py', 'analyze_patients.py', 'pivot_external_table.py']
    stages += ['make_patient_figures.py', 'make_compact_profiles.py',
               'make_permutation_panels.py', 'verify_deliverables.py']
    steps = [run_step(work, s) for s in stages]
    if not args.figures_only:
        subprocess.run([sys.executable, str(ROOT/'code/compare_reproduction.py'),
                        '--reference', str(ROOT), '--reproduced', str(work)], check=True, cwd=ROOT)
    (work/'provenance/run_report.json').write_text(json.dumps({
        'mode': 'figures_from_archived_statistics' if args.figures_only else 'offline_numerical_reanalysis',
        'steps': steps, 'source': str(ROOT), 'output': str(work),
        'numeric_reference': 'Frozen CSV values and NPZ arrays compared independently; publication workbooks are not required.',
        'AC_ICAM': 'Archived all-gene summaries and top-30 calls; raw network acquisition not rerun.'
    }, indent=2)+'\n', encoding='utf-8')
    print(f'Finished. Reproduced outputs: {work}')


if __name__ == '__main__':
    main()
