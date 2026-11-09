#!/usr/bin/env python3
"""Export every worksheet without dropping annotation rows or reformatting numbers."""
from __future__ import annotations
import argparse, csv, datetime, hashlib, io, json, os, tempfile
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parents[1]

def atomic(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(dir=path.parent, prefix='.' + path.name + '.')
    with os.fdopen(fd, 'wb') as f:
        f.write(data); f.flush(); os.fsync(f.fileno())
    os.replace(name, path)

def value(v):
    if v is None: return ''
    if isinstance(v, (datetime.date, datetime.datetime, datetime.time)): return v.isoformat()
    return str(v)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workbook', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    a = ap.parse_args()
    workbook = openpyxl.load_workbook(a.workbook, read_only=True, data_only=False)
    expected = [f'ST_{i}' for i in range(1, 15)]
    if workbook.sheetnames != expected:
        raise ValueError(f'Expected all 14 sheets in order, found {workbook.sheetnames}')
    report = {'workbook': a.workbook.name, 'workbook_sha256': hashlib.sha256(a.workbook.read_bytes()).hexdigest(),
              'export_semantics': 'Complete worksheet grids; annotation/header rows retained; formulas exported as formula text; empty cells remain empty.', 'sheets': []}
    for ws in workbook:
        rows = []
        formula_count = 0
        for row in ws.iter_rows(values_only=True):
            row = list(row)
            while row and row[-1] is None: row.pop()
            formula_count += sum(isinstance(v,str) and v.startswith('=') for v in row)
            rows.append([value(v) for v in row])
        buf = io.StringIO(newline=''); writer = csv.writer(buf, lineterminator='\n'); writer.writerows(rows)
        data = buf.getvalue().encode('utf-8'); path = a.output_dir/(ws.title + '.csv'); atomic(path, data)
        with path.open(newline='', encoding='utf-8') as f:
            if list(csv.reader(f)) != rows: raise AssertionError(f'Round-trip failure: {ws.title}')
        report['sheets'].append({'sheet':ws.title,'file':path.name,'rows':len(rows),
            'max_columns':max(map(len, rows),default=0),'formula_cells':formula_count,
            'sha256':hashlib.sha256(data).hexdigest(),'round_trip_verified':True})
    workbook.close()
    report['passed'] = all(x['round_trip_verified'] for x in report['sheets'])
    atomic(a.output_dir/'export_manifest.json', (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'passed':report['passed'], 'worksheets':len(report['sheets']),
                      'rows':sum(x['rows'] for x in report['sheets'])}))

if __name__ == '__main__': main()
