"""Generate three single-sheet supplemental workbooks using artifact_tool only."""
import os
os.environ['ARTIFACT_TOOL_RPC_DAEMON_STARTUP_TIMEOUT_S']='90'
from pathlib import Path
import csv,re,json
from artifact_tool import Workbook,SpreadsheetFile
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'tables';OUT.mkdir(exist_ok=True)

def col(n):
 s=''
 while n:n,r=divmod(n-1,26);s=chr(65+r)+s
 return s

def convert(v):
 if v=='':return None
 if v in ['True','False']:return 'YES' if v=='True' else 'NO'
 if re.fullmatch(r'-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?',v):
  x=float(v);return int(x) if x.is_integer() and abs(x)<1e15 else x
 return v

CONFIGS=[
 ('TCGA_primary','primary_466_results.csv','Supplementary_Table_7_TCGA_466.xlsx',
 ['DNA_display_rank','Gene','DNA_driver_CNA_n_profiled','DNA_driver_CNA_altered','DNA_driver_CNA_frequency','DNA_driver_CNA_CI95_low','DNA_driver_CNA_CI95_high','DNA_top30_boot_frequency','Functional_pattern','Driver_mutation_frequency','Amplification_frequency','Deep_deletion_frequency'],
 'Patient alterations in 466 COREAD dependency candidates',
 ['Primary: COAD-GDC + READ-GDC | 631 patients; 572 jointly profiled for mutation and CNA (421 colon, 151 rectum).',
 'DNA rank = unique patients with an export-annotated putative-driver mutation OR amplification OR deep deletion; frequencies are not added.',
 'All 466 candidates retained. Profile-based denominators are not verified gene-level callability. A sensitivity excludes the 3 patients represented by 2 samples.',
 'RNA High/Low are source-specific exported calls. No tumor-normal inference or protein result. Missing measurements remain blank.',
 'Rank stability: 2,000 stratified patient-bootstrap draws. Comparative P/FDR describe cohort differences, not candidate-set error control.']),
 ('External_and_robustness','external_466_results_wide.csv','Supplementary_Table_8_External_466.xlsx',
 ['Gene','Primary_DNA_display_rank','Cohort','Cohort_role','Driver_status','DNA_driver_CNA_N','DNA_driver_CNA_Altered','DNA_driver_CNA_Frequency','DNA_driver_CNA_TCGA_frequency','DNA_driver_CNA_Heterogeneity_P','DNA_driver_CNA_Heterogeneity_BH_FDR'],
 'External alteration evidence and paired TCGA robustness',
 ['466 genes x 3 comparison cohorts (1,398 rows). Filter Cohort and Gene; endpoints occupy separate column blocks. AC-ICAM and CPTAC-2 are external colon cohorts; PanCancer reuses TCGA patients.',
 'AC-ICAM: 281 DNA patients (456 exact symbols); 348 RNA patients (451 exact symbols). Driver-specific labels are unavailable; broad mutations are not drivers.',
 'CPTAC-2: 110 patients, 105 jointly DNA-profiled. Driver OR CNA and broad coding/splice OR CNA are reported as separate endpoints.',
 'Fisher P/FDR test rate HETEROGENEITY, not successful validation. PanCancer: paired exact McNemar and matched-patient denominators.',
 'RNA thresholds/references not harmonized: descriptive only. Protein gene-level coverage unavailable. No absent measurement is assigned a zero alteration frequency.']),
 ('Permutation_statistics','permutation_statistics.csv','Supplementary_Table_9_Permutations.xlsx',
 ['Test_ID','Reference','Comparison','Endpoint','Null_scheme','Matched_genes','Units_tested','Spearman_rho','Permutations','Null_exceedances','Empirical_P_greater','BH_FDR','BH_family','Status'],
 'Global concordance: 10,000 permutations per evaluable test',
 ['50 planned contrasts: 46 evaluable tests, each with 10,000 permutations. Four driver-specific AC-ICAM contrasts are unavailable, not negative.',
 'Spearman concordance uses all matched genes. Positive alignment test: P = (b + 1)/(B + 1); never P = 0. These are not per-gene validation tests.',
 'BH families: 18 primary external tests; 18 colon-only external sensitivities; 10 paired-TCGA robustness tests.',
 'CNA sensitivity aggregates matched genes by chromosome arm (at least 3 genes per arm) before permutation to reduce regional gene pseudoreplication.',
 'The ZIP preserves every null draw and the 2,000 patient-bootstrap rank draws. Statistical significance does not prove equivalent frequencies or biological causality.'])]
NOTES=[
 'The 466-gene list was frozen from the functional dependency union. No new alteration-frequency cutoff is used to remove genes.',
 'Mutation labels: explicit (putative driver) plus splice_rec. Explicit-label-only sensitivity is included. No uniform fresh variant annotation was generated.',
 'CNA: amplification and deep deletion are kept separate. Composite DNA counts each patient once and requires both mutation and CNA profile flags.',
 'The primary patient-level exports aggregate multiple samples in 3 patients. Single-sample-only sensitivity is separate; a uniform primary-tumor specimen cannot be verified.',
 'Wilson 95% confidence intervals are pointwise. The rank bootstrap resamples patient columns within COAD and READ, not independent genes.',
 'AC-ICAM CNA and mutation sample sets contain the same 281 unique patients. Exact-symbol matching leaves 10 candidates unavailable. Broad mutation = specified coding/splice MAF classes.',
 'AC-ICAM RNA uses source all-sample z-scores with strict z > 2 and z < -2. Supplied TCGA/CPTAC export thresholds are undocumented, so RNA is not formally validated across studies.',
 'PanCancer part 2 has only 524 patients. Its missing 70 patients are not treated as negative. RNA and RPPA coverage cannot be verified at gene level.',
 'Gene-level Fisher tests assess between-cohort rate differences. BH families include all eligible genes x DNA endpoints for each comparison, before ranking.',
 'Global permutation P uses the number of null statistics at least as large as the observed positive Spearman correlation, with plus-one correction. CNA arm-level tests are sensitivity analyses, not definitive correction for every genomic dependency.',
 'Source URLs, code, sample-level categorical calls, all null draws, checksums, and interpretation limits are preserved in the companion reproducibility ZIP.'
]
checks=[]
for key,fn,dest,front,title,notes in CONFIGS:
 if os.environ.get('TABLE_ONLY') and key != os.environ['TABLE_ONLY']:continue
 print('START',key,flush=True)
 records=list(csv.DictReader((ROOT/'results'/fn).open()))
 headers=front+[h for h in records[0] if h not in front]
 if key=='External_and_robustness':records.sort(key=lambda r:({'CPTAC2':0,'AC_ICAM':1,'PANCAN':2}[r['Cohort']],int(r['Primary_DNA_display_rank'])))
 values=[[convert(r[h]) for h in headers] for r in records];n=len(values);end=8+n;last=col(len(headers));idx={h:col(i+1) for i,h in enumerate(headers)}
 wb=Workbook.create();sh=wb.worksheets.add(key)
 sh.get_range(f'A1:{last}{end}').format={'font':{'size':10},'row_height':21,'column_width':16,'vertical_alignment':'center'}
 sh.merge_cells('A1:L1');sh.get_range('A1').values=[['Supplementary Table | '+title]]
 sh.get_range('A1:L1').format={'fill':'#17365D','font':{'bold':True,'color':'#FFFFFF','size':15},'row_height':34}
 for r,t in enumerate(notes,2):
  sh.merge_cells(f'A{r}:L{r}');sh.get_range(f'A{r}').values=[[t]];sh.get_range(f'A{r}:L{r}').format={'fill':'#E9EFF6','wrap_text':True,'row_height':29}
 sh.get_range_by_indexes(7,0,1,len(headers)).values=[headers]
 sh.get_range(f'A8:{last}8').format={'fill':'#24476A','font':{'bold':True,'color':'#FFFFFF','size':10},'row_height':58,'wrap_text':True}
 for a in range(0,n,2000):sh.get_range_by_indexes(8+a,0,min(2000,n-a),len(headers)).values=values[a:a+2000]
 formulas=[]
 if key=='TCGA_primary':
  for ep in ['DNA_driver_CNA','Driver_mutation','Amplification','Deep_deletion','CNA_any','Any_exported_mutation','DNA_broad','RNA_high','RNA_low','Protein_high','Protein_low']:
   formulas.append((ep+'_frequency',ep+'_n_profiled',ep+'_altered'))
  for c in ['COAD_GDC','READ_GDC']:formulas.append((c+'_DNA_frequency',c+'_DNA_n',c+'_DNA_altered'))
 elif key=='External_and_robustness':formulas=[(ep+'_Frequency',ep+'_N',ep+'_Altered') for ep in ['DNA_driver_CNA','Driver_mutation','Amplification','Deep_deletion','CNA_any','Any_exported_mutation','DNA_broad','RNA_high','RNA_low']]
 for d,nf,kf in formulas:
  sh.get_range(f'{idx[d]}9:{idx[d]}{end}').formulas=[[f'=IF({idx[nf]}{r}>0,{idx[kf]}{r}/{idx[nf]}{r},"")'] for r in range(9,end+1)]
 if key=='Permutation_statistics':sh.get_range(f'{idx["Empirical_P_greater"]}9:{idx["Empirical_P_greater"]}{end}').formulas=[[f'=IF({idx["Permutations"]}{r}>0,({idx["Null_exceedances"]}{r}+1)/({idx["Permutations"]}{r}+1),"")'] for r in range(9,end+1)]
 for h,c in idx.items():
  rg=sh.get_range(f'{c}9:{c}{end}')
  if ('frequency' in h.lower() or 'CI95' in h) and not any(s in h for s in ['rank','Chronos','tail']):rg.set_number_format('0.0%')
  elif h.endswith('_P') or h.endswith('_FDR') or h in ['Empirical_P_greater','BH_FDR']:rg.set_number_format('0.00E+00')
  elif h in ['Spearman_rho','Null_mean','Null_SD','Organoid_Chronos_mean']:rg.set_number_format('0.000')
  elif 'rank' in h.lower():rg.set_number_format('0.0')
  if h=='Gene':sh.get_range(f'{c}:{c}').format.column_width=14
  if any(t in h for t in ['Interpretation','qualification','Source_URL','Null_scheme','Cohort_role','Status','BH_family']):
   sh.get_range(f'{c}:{c}').format.column_width=36;rg.format.wrap_text=True
 sh.tables.add(f'A8:{last}{end}',True,key+'_table');sh.freeze_panes.freeze_rows(8);sh.freeze_panes.freeze_columns(2)
 nr=end+3;sh.merge_cells(f'A{nr}:L{nr}');sh.get_range(f'A{nr}').values=[['METHODS AND INTERPRETATION']];sh.get_range(f'A{nr}:L{nr}').format={'fill':'#17365D','font':{'bold':True,'color':'#FFFFFF'},'row_height':25}
 for r,t in enumerate(NOTES,nr+1):sh.merge_cells(f'A{r}:L{r}');sh.get_range(f'A{r}').values=[[t]];sh.get_range(f'A{r}:L{r}').format={'wrap_text':True,'row_height':30}
 print('VERIFY',key,flush=True)
 ck=wb.inspect({'kind':'table','range':f'{key}!A8:L12','include':'values,formulas','table_max_rows':5,'table_max_cols':12})
 (ROOT/'provenance'/f'{key}_inspection.txt').write_text(ck.ndjson)
 err=wb.inspect({'kind':'match','search_term':'#REF!|#DIV/0!|#VALUE!|#NAME\\?','options':{'use_regex':True,'max_results':20}})
 (ROOT/'provenance'/f'{key}_formula_scan.txt').write_text(err.ndjson)
 print('RENDER',key,flush=True)
 try:wb.render({'sheet_name':key,'range':'A1:L15','scale':1}).save(str(ROOT/'provenance'/f'{key}_preview.png'))
 except Exception as ex:print('RENDER WARNING',type(ex).__name__,str(ex),flush=True)
 print('EXPORT',key,flush=True)
 SpreadsheetFile.export_xlsx(wb).save(str(OUT/dest))
 checks.append({'file':dest,'sheets':1,'rows':n,'columns':len(headers),'source_csv':fn})
 print('SAVED',dest,flush=True)
(ROOT/'provenance'/('workbook_build_'+os.environ.get('TABLE_ONLY','all')+'.json')).write_text(json.dumps(checks,indent=2))
print('ALL WORKBOOKS SAVED',flush=True)
