"""Compare a clean rerun with the shipped numerical reference, without Excel engines."""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np

CSV_FILES = ['primary_466_results.csv', 'external_466_results_long.csv',
             'external_466_results_wide.csv', 'permutation_statistics.csv',
             'colon_rectum_heterogeneity.csv', 'exported_events_long.csv']
NPZ_FILES = ['COAD_GDC_calls.npz', 'READ_GDC_calls.npz', 'COREAD_GDC_calls.npz',
             'CPTAC2_calls.npz', 'PANCAN_P1_calls.npz', 'PANCAN_P2_calls.npz',
             'PANCAN_calls.npz', 'primary_rank_bootstrap.npz', 'permutation_null_draws.npz']


def compare(reference: Path, reproduced: Path) -> dict:
    checks = []
    for name in CSV_FILES:
        with (reference/'results'/name).open(newline='', encoding='utf-8') as f:
            a = list(csv.reader(f))
        with (reproduced/'results'/name).open(newline='', encoding='utf-8') as f:
            b = list(csv.reader(f))
        if len(a) != len(b) or a[0] != b[0]:
            raise AssertionError(f'{name}: row count or headers differ')
        count = 0
        for i, (ra, rb) in enumerate(zip(a[1:], b[1:]), 2):
            if len(ra) != len(rb):
                raise AssertionError(f'{name}: row {i} widths differ')
            for j, (x, y) in enumerate(zip(ra, rb)):
                if x == y:
                    count += 1
                    continue
                try:
                    equal = math.isclose(float(x), float(y), abs_tol=1e-12, rel_tol=1e-11)
                except (TypeError, ValueError):
                    equal = False
                if not equal:
                    raise AssertionError(f'{name} row {i} {a[0][j]}: {x!r} != {y!r}')
                count += 1
        checks.append({'file': name, 'rows': len(a)-1, 'cells_compared': count, 'discrepancies': 0})
    for name in NPZ_FILES:
        a = np.load(reference/'results'/name, allow_pickle=False)
        b = np.load(reproduced/'results'/name, allow_pickle=False)
        if a.files != b.files:
            raise AssertionError(f'{name}: named arrays differ')
        maxdiff = 0.
        for k in a.files:
            x, y = a[k], b[k]
            if x.shape != y.shape:
                raise AssertionError(f'{name}:{k} shape mismatch')
            if x.dtype.kind in 'fciub':
                if not np.allclose(x, y, rtol=1e-11, atol=1e-12, equal_nan=True):
                    raise AssertionError(f'{name}:{k} values differ')
                if x.size and x.dtype.kind in 'fc':
                    maxdiff = max(maxdiff, float(np.nanmax(np.abs(x-y))))
            elif not np.array_equal(x, y):
                raise AssertionError(f'{name}:{k} identities differ')
        checks.append({'file': name, 'arrays_compared': len(a.files),
                       'max_abs_float_difference': maxdiff, 'discrepancies': 0})
    report = {'result': 'PASS', 'rtol': 1e-11, 'atol': 1e-12,
              'checks': checks, 'scope': 'CSV values and NPZ arrays; does not verify raw sequencing or biological validity.'}
    (reproduced/'provenance/numerical_regression.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--reproduced', type=Path, required=True)
    args = p.parse_args()
    r = compare(args.reference, args.reproduced)
    print('PASS: clean rerun agrees with reference for', len(r['checks']), 'numerical files')
