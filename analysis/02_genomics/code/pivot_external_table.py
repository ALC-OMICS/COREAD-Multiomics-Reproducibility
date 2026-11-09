from pathlib import Path
from collections import defaultdict
import csv
R=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((R/'results/external_466_results_long.csv').open()));groups=defaultdict(dict)
for r in rows:groups[(r['Gene'],r['Cohort'])][r['Endpoint']]=r
out=[]
dna=['DNA_driver_CNA','Driver_mutation','Amplification','Deep_deletion','CNA_any','Any_exported_mutation','DNA_broad']
for (g,c),rs in groups.items():
 first=rs['DNA_driver_CNA'];o={'Gene':g,'Primary_DNA_display_rank':first['Primary_DNA_display_rank'],'Cohort':c,'Cohort_role':first['Cohort_role'],'Driver_status':rs['Driver_mutation']['Status']}
 for ep in dna:
  r=rs[ep]
  for dest,src in [('N','External_n_evaluable'),('Altered','External_altered'),('Frequency','External_frequency'),('CI95_low','External_CI95_low'),('CI95_high','External_CI95_high'),('TCGA_frequency','TCGA_frequency'),('Heterogeneity_P','Heterogeneity_P'),('Heterogeneity_BH_FDR','Heterogeneity_BH_FDR'),('COAD_heterogeneity_P','COAD_heterogeneity_P'),('COAD_heterogeneity_BH_FDR','COAD_heterogeneity_BH_FDR')]:o[ep+'_'+dest]=r[src]
 for ep in ['RNA_high','RNA_low']:
  r=rs[ep]
  for dest,src in [('N','External_n_evaluable'),('Altered','External_altered'),('Frequency','External_frequency'),('CI95_low','External_CI95_low'),('CI95_high','External_CI95_high')]:o[ep+'_'+dest]=r[src]
 o['RNA_comparison_status']='SOURCE_SPECIFIC_NOT_HARMONIZED' if c!='PANCAN' else 'NO_VERIFIED_QUERY_PROFILE_DENOMINATOR'
 o['Protein_high_observed_calls']=rs['Protein_high']['Observed_export_calls_no_coverage_claim'] if c=='PANCAN' else ''
 o['Protein_low_observed_calls']=rs['Protein_low']['Observed_export_calls_no_coverage_claim'] if c=='PANCAN' else ''
 o['Protein_status']='NO_VERIFIED_GENE_PROFILE_COVERAGE'
 for k in ['Paired_n','Paired_TCGA_altered','Paired_PanCancer_altered','TCGA_only_call_n','PanCancer_only_call_n','Paired_call_agreement']:o['DNA_'+k]=first[k]
 o['Driver_interpretation']='Export putative-driver labels; annotation source/version not independently verified.' if c!='AC_ICAM' else 'No driver labels in retrieved MAF; broad coding/splice counts are NOT drivers.'
 o['Source_URL']=first['Source_URL'];out.append(o)
out.sort(key=lambda r:({'CPTAC2':0,'AC_ICAM':1,'PANCAN':2}[r['Cohort']],int(r['Primary_DNA_display_rank'])))
with (R/'results/external_466_results_wide.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
print(len(out),len(out[0]))
