#!/usr/bin/env python3
"""Combine the 31 panel maps and make lossless, auditable figure-level TSV containers."""
from __future__ import annotations
import argparse,csv,hashlib,io,json,os,re,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
FIELDS=['figure','panel','source_files','description','units','missingness']
def atomic(p,data):
    p.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(dir=p.parent,prefix='.'+p.name+'.')
    with os.fdopen(fd,'wb')as f:f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(tmp,p)
def rows(path):
    with path.open(newline='',encoding='utf-8-sig')as f:
        return list(csv.reader(f,delimiter='\t'if path.suffix.lower()=='.tsv'else','))
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,default=ROOT/'source_data');a=ap.parse_args()
    panels=[]
    for p in sorted((ROOT/'source_data').glob('*/panel_index.csv')):
        with p.open(newline='',encoding='utf-8-sig')as f:
            for row in csv.DictReader(f):
                row={k:row[k]for k in FIELDS};number=int(re.search(r'\d+',row['figure']).group())
                row['figure']=f'Figure {number}';panels.append(row)
    panels.sort(key=lambda r:(int(r['figure'].split()[-1]),r['panel']))
    expected={1:'ABCDE',2:'ABCD',3:'ABCD',4:'ABCD',5:'ABCD',6:'ABCDEF',7:'ABCD'}
    if len(panels)!=31:raise AssertionError(f'Expected 31 panels, found {len(panels)}')
    for n,letters in expected.items():
        if ''.join(x['panel']for x in panels if x['figure']==f'Figure {n}')!=letters:raise AssertionError(f'Panel map differs for Figure {n}')
    b=io.StringIO(newline='');w=csv.DictWriter(b,fieldnames=FIELDS,lineterminator='\n');w.writeheader();w.writerows(panels)
    atomic(a.output_dir/'panel_index.csv',b.getvalue().encode())
    report={'format':'COREAD figure-level TSV container, version 1','figures':[],'panels':31,'passed':False}
    for n in expected:
        records=[['FORMAT','COREAD_SOURCE_DATA','1'],['FIGURE',str(n)]];truth=[];sources=[]
        for panel in (p for p in panels if p['figure']==f'Figure {n}'):
            records.append(['PANEL',panel['panel'],panel['description'],panel['units'],panel['missingness']])
            for rel in panel['source_files'].split(';'):
                rel=rel.strip();p=ROOT/rel
                if not p.is_file():raise FileNotFoundError(p)
                table=rows(p);digest=hashlib.sha256(p.read_bytes()).hexdigest()
                records.append(['TABLE',rel,digest,str(len(table))])
                for i,row in enumerate(table):records.append(['HEADER'if i==0 else'DATA',*row])
                records.append(['END_TABLE']);truth.append(table);sources.append({'path':rel,'sha256':digest,'rows':len(table)})
            records.append(['END_PANEL'])
        out=io.StringIO(newline='');csv.writer(out,delimiter='\t',lineterminator='\n').writerows(records)
        data=out.getvalue().encode('utf-8');path=a.output_dir/f'SourceData_Figure_{n}.tsv';atomic(path,data)
        decoded=[];current=None
        with path.open(newline='',encoding='utf-8')as f:
            for rec in csv.reader(f,delimiter='\t'):
                if rec[0]=='TABLE':current=[]
                elif rec[0]in('HEADER','DATA'):current.append(rec[1:])
                elif rec[0]=='END_TABLE':decoded.append(current);current=None
        if decoded!=truth:raise AssertionError(f'Lossless round-trip failed for Figure {n}')
        report['figures'].append({'figure':n,'file':path.name,'sha256':hashlib.sha256(data).hexdigest(),
            'tables':sources,'total_table_rows':sum(map(len,truth)),'round_trip_verified':True})
    report['passed']=True;atomic(a.output_dir/'source_data_manifest.json',(json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps({'passed':True,'figures':7,'panels':31,'table_blocks':sum(len(f['tables'])for f in report['figures'])}))
if __name__=='__main__':main()
