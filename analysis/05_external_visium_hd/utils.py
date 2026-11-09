from pathlib import Path
import json, os
import numpy as np
import pandas as pd

def csv_write(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('w', encoding='utf-8', newline='') as f:
        frame.to_csv(f, index=False)
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)


def json_write(obj, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp')
    with temp.open('w', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ix = np.flatnonzero(np.isfinite(p))
    order = ix[np.argsort(p[ix], kind='stable')]
    if len(order):
        adj = p[order] * len(order) / np.arange(1, len(order) + 1)
        q[order] = np.minimum(1, np.minimum.accumulate(adj[::-1])[::-1])
    return q


def flag(values):
    return values.astype(str).str.lower().isin(['true', 'yes', '1'])


def compare_frames(actual, expected, keys=('Gene',), rtol=1e-9, atol=1e-11):
    a = actual.set_index(list(keys)).sort_index()
    e = expected.set_index(list(keys)).sort_index()
    failures = []
    if not a.index.equals(e.index):
        return {'passed': False, 'failures': ['key mismatch'], 'checked_columns': 0}
    for col in e:
        if col not in a:
            failures.append({'column': col, 'reason': 'missing'})
            continue
        ev, av = e[col], a[col]
        if pd.api.types.is_bool_dtype(ev):
            equal = flag(av).to_numpy() == flag(ev).to_numpy()
        elif pd.api.types.is_numeric_dtype(ev):
            x = pd.to_numeric(av, errors='coerce').to_numpy(float)
            y = ev.to_numpy(float)
            equal = np.isclose(x, y, rtol=rtol, atol=atol, equal_nan=True)
        else:
            equal = av.fillna('').astype(str).to_numpy() == ev.fillna('').astype(str).to_numpy()
        if not np.all(equal):
            bad = np.flatnonzero(~equal)
            failures.append({'column': col, 'n_mismatch': len(bad),
                             'example_keys': [str(e.index[i]) for i in bad[:5]]})
    return {'passed': not failures, 'checked_columns': len(e.columns), 'rows': len(e), 'failures': failures}

