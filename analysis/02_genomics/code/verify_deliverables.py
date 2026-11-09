"""Verify numerical invariants without requiring editorial Excel workbooks."""
from pathlib import Path
import csv,json,math
import numpy as np
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parents[1]
report=[]
for name,n in [('primary_466_results.csv',466),('external_466_results_wide.csv',1398),('permutation_statistics.csv',50)]:
    with (R/'results'/name).open(newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==n,(name,len(rows),n)
    keys=[(r.get('Cohort',''),r.get('Gene',r.get('Test_ID'))) for r in rows]
    assert len(set(keys))==n,'Duplicate numerical record: '+name
    report.append({'file':name,'records':n,'unique_identifiers':True})
stat=list(csv.DictReader((R/'results/permutation_statistics.csv').open()));null=np.load(R/'results/permutation_null_draws.npz');fams={}
for t in stat:
 if t['Status']!='EVALUATED_CONDITIONAL_CONCORDANCE':continue
 a=null[t['Test_ID']];assert len(a)==10000;b=int(np.sum(a>=float(t['Spearman_rho'])-1e-12));assert b==int(t['Null_exceedances']);assert math.isclose((b+1)/10001,float(t['Empirical_P_greater']),rel_tol=1e-12)
 fams.setdefault(t['BH_family'],[]).append(t)
for family,ts in fams.items():
 q=false_discovery_control([float(t['Empirical_P_greater']) for t in ts],method='bh')
 assert np.allclose(q,[float(t['BH_FDR']) for t in ts],rtol=1e-12,atol=1e-14)
for p in (R/'figures').glob('*.svg'):
 text=p.read_text();assert '<image' not in text,'Raster in '+p.name;assert '<text' in text,'No editable text '+p.name
(R/'provenance/final_verification.json').write_text(json.dumps({'numerical_records':report,'permutation_tests_verified':len(null.files),'publication_workbooks_required':False,'null_draws_each':10000,'BH_independent_check':'scipy.stats.false_discovery_control','svg_figures_verified':len(list((R/'figures').glob('*.svg')))},indent=2))
print(json.dumps(report,indent=2));print('VERIFIED',len(null.files),'permutation tests and all SVG figures')
