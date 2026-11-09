"""Verify repository integrity and numerical agreement with the frozen numerical CSV references."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from common import ROOT,load_csv,compare_rows,number

def run(build=None, hashes=True):
    report={}
    manifest=ROOT/'provenance/sha256_manifest.json'
    if hashes and manifest.exists():
        items=json.loads(manifest.read_text());bad=[]
        for path,digest in items.items():
            p=ROOT/path
            if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:bad.append(path)
        if bad:raise AssertionError('Checksum mismatch: '+str(bad[:10]))
        report['checksums_verified']=len(items)
    references=['cell_line_results_1165.csv','integration_results_1165.csv','organoid_results_1720.csv']
    for fname in references:
        expected=load_csv(ROOT/'results/csv'/fname)
        if len(expected)!=1720: raise AssertionError('Incomplete numeric reference: '+fname)
        report[fname]={'reference_rows':len(expected),'reference_columns':len(expected[0])}
        if build:
            sub='organoids' if fname.startswith('organoid_') else 'cell_lines'
            calculated=load_csv(Path(build)/sub/fname)
            report[fname]['recalculated_cells_checked']=compare_rows(calculated,expected)
    cl=load_csv(ROOT/'results/csv/cell_line_results_1165.csv');org=load_csv(ROOT/'results/csv/organoid_results_1720.csv');it=load_csv(ROOT/'results/csv/integration_results_1165.csv')
    genes=[r['Gene_symbol'] for r in cl];assert len(genes)==len(set(genes))==1720
    assert genes==[r['Gene_symbol'] for r in org]==[r['Gene_symbol'] for r in it]
    counts={'cell_line_candidates':sum(r['CL_candidate']=='YES' for r in cl),'organoid_candidates':sum(r['Primary_candidate_3DO']=='YES' for r in org),'union':sum(r['Union_candidate']=='YES' for r in it),'shared':sum(r['Both_systems_candidate']=='YES' for r in it),'both_evaluable':sum(r['Both_evaluable']=='YES' for r in it),'higher_frequency':sum(r['Frequency_preference']=='HIGHER_IN_COREAD' for r in cl)}
    assert counts=={'cell_line_candidates':335,'organoid_candidates':324,'union':466,'shared':193,'both_evaluable':1613,'higher_frequency':23}
    report['counts']=counts
    report['status']='PASS'
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--build',type=Path);p.add_argument('--no-hashes',action='store_true');a=p.parse_args()
    report=run(a.build,not a.no_hashes);print(json.dumps(report,indent=2))
    if a.build:(a.build/'verification_report.json').write_text(json.dumps(report,indent=2))
