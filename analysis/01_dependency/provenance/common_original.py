"""Small I/O utilities used by the reproducibility entrypoints."""
from pathlib import Path
import csv, json, re, math
import numpy as np
from xlsx_snapshot import parse_rows
ROOT = Path(__file__).resolve().parents[1]

def load_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def read_gene_table(path):
    """Read the contiguous input table; no inference from note/metadata rows."""
    raw = {i: {re.match('[A-Z]+', k)[0]: v for k, v in row.items()} for i, row in parse_rows(path, 1)}
    head = next(i for i, r in raw.items() if r.get('A') == 'Input_order')
    headers = list(raw[head].values()); columns = list(raw[head])
    records = []
    for i in range(head + 1, head + 1721):
        records.append({h: raw[i].get(c) for c,h in zip(columns,headers)})
    assert [int(r['Input_order']) for r in records] == list(range(1,1721))
    return headers, records, raw

def number(x):
    if x is None or x == '': return np.nan
    return float(x)

def clean(x):
    if isinstance(x,np.generic): x=x.item()
    if isinstance(x,float):
        if np.isnan(x): return None
        if np.isposinf(x): return 'Inf'
        if np.isneginf(x): return '-Inf'
    return x

def write_csv(path, headers, records):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=headers);w.writeheader()
        for r in records:w.writerow({h:clean(r.get(h)) for h in headers})

def bh(p):
    p=np.asarray(p,float); q=np.full_like(p,np.nan); ok=np.isfinite(p)
    values=p[ok]; order=np.argsort(values,kind='stable'); m=len(values)
    if m:
        adj=np.minimum.accumulate((values[order]*m/np.arange(1,m+1))[::-1])[::-1]
        result=np.empty(m);result[order]=np.minimum(adj,1);q[ok]=result
    return q

def compare_rows(actual, expected, columns=None, rtol=1e-8, atol=1e-11):
    if len(actual)!=len(expected): raise AssertionError('Different number of rows')
    differences=[];checked=0
    for idx,(a,e) in enumerate(zip(actual,expected)):
        for k in (columns or e.keys()):
            checked+=1; x=a.get(k);y=e.get(k)
            if (x is None or x=='') and (y is None or y==''):continue
            if x==y:continue
            try:
                nx, ny = number(x), number(y)
                if (math.isnan(nx) and math.isnan(ny)) or nx == ny:
                    continue
                if math.isfinite(nx) and math.isfinite(ny) and abs(nx - ny) <= atol + rtol * abs(ny):
                    continue
            except (ValueError,TypeError):pass
            differences.append({'row':idx+1,'column':k,'actual':str(x),'expected':str(y)})
    if differences:raise AssertionError(json.dumps(differences[:12],indent=2))
    return checked
