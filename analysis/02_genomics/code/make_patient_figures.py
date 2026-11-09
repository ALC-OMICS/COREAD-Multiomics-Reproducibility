"""Exact-data patient figures. Each figure contains one plot, with editable SVG text.
Default matplotlib colour cycle/colormap only. No images or reconstructed values.
"""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.ticker import PercentFormatter
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams['svg.fonttype']='none';plt.rcParams['font.family']='DejaVu Sans';plt.rcParams['font.size']=10
summary=json.loads((ROOT/'results/analysis_summary.json').read_text());S=json.loads((ROOT/'results/all_endpoint_summaries.json').read_text())
genes=(ROOT/'input/union_candidates_466.txt').read_text().split();gi={g:i for i,g in enumerate(genes)};top=summary['primary_top30'];ti=np.array([gi[g] for g in top])
stats=list(csv.DictReader((ROOT/'results/permutation_statistics.csv').open()));primary=list(csv.DictReader((ROOT/'results/primary_466_results.csv').open()));legends=[]

def finalize(fig,stem,legend):
 fig.savefig(OUT/(stem+'.svg'),bbox_inches='tight',pad_inches=.12)
 fig.savefig(OUT/(stem+'.png'),dpi=180,bbox_inches='tight',pad_inches=.12)
 plt.close(fig);legends.append(stem+'\n'+legend+'\n');print('SAVED',stem,flush=True)

def oncoprint(cohort,stem,ac=False):
 if ac:
  z=np.load(ROOT/'results/AC_ICAM_top30_calls.npz');assert list(z['genes'])==top
  codes=z['codes'];patients=z['patients'];N=len(patients)
  events=[(codes&v)>0 for v in [1,2,4,8,16]];rna_avail=(codes&32)==0
  dna=np.any(np.stack(events[:3]),axis=0);counts=dna.sum(1);strata=[np.arange(N)]
  mutation_name='Coding/splice\nmutation*';title='AC-ICAM | Source-specific external OncoPrint'
  subtitle='Same 30 genes selected in TCGA; 281 patients. *Mutation calls are NOT driver-annotated.'
 else:
  z=np.load(ROOT/f'results/{cohort}_calls.npz');joint=z['avail__MUTATIONS']&z['avail__CNA'];ix=np.flatnonzero(joint);patients=z['patients'][ix];N=len(ix)
  events=[z['event__'+e][ti][:,ix].astype(bool) for e in ['Driver_mutation','Amplification','Deep_deletion','RNA_high','RNA_low']]
  rna_avail=np.broadcast_to(z['avail__MRNA'][ix],(30,N));dna=np.any(np.stack(events[:3]),axis=0);counts=dna.sum(1)
  strata=[np.flatnonzero(ix<461),np.flatnonzero(ix>=461)] if cohort=='COREAD_GDC' else [np.arange(N)]
  mutation_name='Putative-driver\nmutation';title='TCGA COREAD | Patient-level alteration landscape' if cohort=='COREAD_GDC' else 'CPTAC-2 | External colon cohort'
  subtitle=f'30 genes fixed by primary DNA rank; {N} patients jointly profiled for mutation and CNA.'
 def ordergroup(ids):
  return sorted(ids,key=lambda j:tuple(-int(x) for x in dna[:,j])+(-int(dna[:,j].sum()),str(patients[j])))
 order=np.array([j for ids in strata for j in ordergroup(ids)])
 events=[e[:,order] for e in events];rna_avail=rna_avail[:,order];patients=patients[order]
 fig=plt.figure(figsize=(12.8,9.6));ax=fig.add_axes([.12,.18,.75,.70])
 patches=[];values=[]
 # Five narrow lanes per gene preserve co-occurring alterations without stacking frequencies.
 for k,e in enumerate(events):
  for i,j in np.argwhere(e):patches.append(Rectangle((j,i+0.10+0.15*k),1,.14));values.append(k+1)
 pc=PatchCollection(patches,linewidths=0);pc.set_array(np.array(values));pc.set_clim(1,5);ax.add_collection(pc)
 for i,j in np.argwhere(~rna_avail):ax.text(j+.5,i+.76,'×',fontsize=2.2,ha='center',va='center')
 ax.set_xlim(0,N);ax.set_ylim(30,0);ax.set_yticks(np.arange(30)+.45);ax.set_yticklabels(top,fontstyle='italic',fontsize=8.3);ax.tick_params(axis='y',length=0,pad=5)
 ax.set_xticks([]);ax.spines[['left','top','right']].set_visible(False)
 if len(strata)==2:
  boundary=len(strata[0]);ax.axvline(boundary,linewidth=.6,linestyle=':')
  ax.text(boundary/2,1.015,f'Colon (n = {boundary})',transform=ax.get_xaxis_transform(),ha='center',fontsize=9)
  ax.text((N+boundary)/2,1.015,f'Rectum (n = {N-boundary})',transform=ax.get_xaxis_transform(),ha='center',fontsize=9)
 for i,k in enumerate(counts):ax.text(1.012,i+.45,f'{100*k/N:.1f}% ({k}/{N})',transform=ax.get_yaxis_transform(),fontsize=7.7,va='center')
 ax.text(1.012,-.35,'Broad DNA*' if ac else 'DNA union',transform=ax.get_yaxis_transform(),fontsize=8,fontweight='bold')
 ax.set_xlabel('Patients; all joint-profile cases retained, including cases without a displayed DNA alteration',fontsize=9,labelpad=10)
 cb=fig.colorbar(pc,ax=ax,orientation='horizontal',fraction=.035,pad=.08,aspect=55,boundaries=np.arange(.5,6,1),values=np.arange(1,6));cb.solids.set_rasterized(False);cb.set_ticks([1,2,3,4,5]);cb.set_ticklabels([mutation_name,'Amplification','Deep deletion','mRNA high','mRNA low']);cb.ax.tick_params(labelsize=8,length=0)
 fig.text(.12,.96,title,fontsize=14,fontweight='bold');fig.text(.12,.931,subtitle,fontsize=9.5)
 fig.text(.12,.065,'Each gene has five alteration lanes; overlapping events in a patient remain visible. Blank = no exported call in an available profile.',fontsize=8)
 if ac:foot='RNA: source z > 2 / z < −2; reference is cohort-specific. No driver-specific conclusion is supported by this panel.'
 else:foot='RNA: original export High/Low calls, with undocumented thresholds. Unprofiled RNA is marked ×. No protein assay is shown.'
 fig.text(.12,.043,foot,fontsize=8)
 finalize(fig,stem,f'{title}. Gene order is fixed from TCGA DNA union frequencies, not selected in the external cohort. {subtitle} '+foot+' Absence of a call is conditional on source profile flags; gene callability was not independently verified.')
 codesout=sum((events[k].astype(np.uint8)*(1<<k) for k in range(5)))+32*(~rna_avail).astype(np.uint8)
 with (OUT/(stem+'_source.csv')).open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['Gene']+list(patients));w.writerows([[g]+list(codesout[i]) for i,g in enumerate(top)])

oncoprint('COREAD_GDC','Figure_1_TCGA_OncoPrint_top30')
oncoprint('CPTAC2','Figure_2_CPTAC2_OncoPrint_top30')
oncoprint('AC_ICAM','Figure_3_ACICAM_OncoPrint_top30',True)

# Per-modality frequencies on separate offsets: no arithmetic sum or normalization.
fig=plt.figure(figsize=(10.8,10.4));ax=fig.add_axes([.12,.13,.72,.75]);y=np.arange(30)
series=[('DNA_driver_CNA','DNA union (joint profiles)'),('Driver_mutation','Putative-driver mutation'),('Amplification','Amplification'),('Deep_deletion','Deep deletion'),('RNA_high','mRNA high (source calls)'),('RNA_low','mRNA low (source calls)')]
for j,(ep,label) in enumerate(series):
 a=S['COREAD_GDC'][ep];f=np.array(a['f'])[ti];ax.scatter(100*f,y+(j-2.5)*.115,s=15,marker=['o','s','^','v','D','P'][j],label=label)
ax.set_ylim(29.65,-.65);ax.set_yticks(y);ax.set_yticklabels(top,fontstyle='italic',fontsize=8.5);ax.set_xlim(-1,100);ax.set_xticks(np.arange(0,101,20));ax.grid(axis='x',linewidth=.4,alpha=.25);ax.set_axisbelow(True);ax.spines[['right','top']].set_visible(False);ax.tick_params(axis='y',length=0)
ax.set_xlabel('Patients with the indicated alteration (%)',fontsize=10)
fig.text(.12,.962,'TCGA COREAD | Frequencies by alteration type',fontsize=14,fontweight='bold')
fig.text(.12,.934,'Same 30 genes as the OncoPrint; no summed marginal frequencies and no unit-sum normalization.',fontsize=9.3)
ax.legend(loc='upper left',bbox_to_anchor=(.995,1.01),frameon=False,fontsize=8,handletextpad=.3)
fig.text(.12,.055,'Denominators: DNA union 572; mutation 581; CNA 621; RNA 624. All denominators are source profile-based.',fontsize=8.5)
fig.text(.12,.033,'RNA High/Low are not tumor-versus-normal differential-expression calls. The DNA ranking does not rank therapeutic actionability.',fontsize=8)
finalize(fig,'Figure_4_TCGA_Frequencies_by_modality','Each point shows a modality-specific patient frequency; categories can overlap and are not summed. Gene order is fixed by the primary DNA union. '+ 'DNA union n=572, mutation n=581, CNA n=621, RNA n=624; source-specific RNA calls are descriptive.')

# Global external rank concordance, all eligible genes; no top-30 filtering.
sel=[r for r in stats if r['Reference']=='COREAD_GDC' and r['Comparison']!='PANCAN' and r['Null_scheme']=='GENE_LABEL_PERMUTATION' and r['Status'].startswith('EVALUATED')]
labels={'DNA_driver_CNA':'Driver / CNA union','Driver_mutation':'Putative-driver mutation','Amplification':'Amplification','Deep_deletion':'Deep deletion','CNA_any':'Amplification / deletion union','Any_exported_mutation':'Broad coding/splice mutations','DNA_broad':'Broad mutation / CNA union'}
fig=plt.figure(figsize=(12,7.0));ax=fig.add_axes([.36,.22,.40,.63]);yy=np.arange(len(sel))
for cohort in ['AC_ICAM','CPTAC2']:
 ix=[j for j,r in enumerate(sel) if r['Comparison']==cohort]
 ax.barh(np.array(ix),[float(sel[j]['Spearman_rho']) for j in ix],height=.60,label=cohort.replace('_','-'))
ax.set_yticks(yy);ax.set_yticklabels([r['Comparison'].replace('_','-')+' | '+labels[r['Endpoint']] for r in sel],fontsize=8.8);ax.invert_yaxis();ax.set_xlim(0,1);ax.set_xlabel('Spearman rank concordance',fontsize=10);ax.set_xticks(np.linspace(0,1,6));ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0);ax.grid(axis='x',linewidth=.4,alpha=.25);ax.set_axisbelow(True)
for j,r in enumerate(sel):
 v=float(r['Spearman_rho']);ax.text(v+.014,j,f'{v:.3f}',va='center',fontsize=8.3)
 ax.text(1.16,j,f"{int(r['Matched_genes'])}   {float(r['Empirical_P_greater']):.2e}   {float(r['BH_FDR']):.2e}",transform=ax.get_yaxis_transform(),va='center',fontsize=8.1)
ax.text(1.16,-.85,'Genes    Empirical P    BH FDR',transform=ax.get_yaxis_transform(),fontsize=8.2,fontweight='bold')
fig.text(.08,.955,'External concordance depends on alteration type',fontsize=14,fontweight='bold')
fig.text(.08,.912,'10,000 gene-label permutations per test; all matched genes. Positive concordance is not equivalent prevalence.',fontsize=9.3)
fig.text(.08,.10,'AC-ICAM driver-specific results are unavailable; broad coding/splice counts are a secondary proxy. RNA is not included in these tests.',fontsize=8.3)
fig.text(.08,.07,'CNA arm-level sensitivities and colon-only comparisons are in the statistical table. PanCancer is separate, non-independent robustness.',fontsize=8.3)
finalize(fig,'Figure_5_External_rank_concordance','Spearman correlations compare all matched genes, not the displayed top 30. Gene-label empirical P values use 10,000 permutations and plus-one correction. BH covers 18 evaluable primary-external tests including CNA arm-level sensitivities. Source annotation and assay differences remain; low effect-size correlations should not be described as strong replication.')

# Permutation distribution for the primary composite in CPTAC: deliberately not
# selecting the strongest available endpoint to illustrate the method.
test=next(r for r in stats if r['Comparison']=='CPTAC2' and r['Reference']=='COREAD_GDC' and r['Endpoint']=='DNA_driver_CNA' and r['Null_scheme']=='GENE_LABEL_PERMUTATION')
null=np.load(ROOT/'results/permutation_null_draws.npz')[test['Test_ID']];obs=float(test['Spearman_rho'])
fig=plt.figure(figsize=(8.5,5.3));ax=fig.add_axes([.12,.20,.81,.65]);ax.hist(null,bins=50,label='10,000 permuted gene correspondences');ax.axvline(obs,linestyle='--',linewidth=1.5,label=f'Observed rho = {obs:.3f}')
ax.set_xlim(min(-.20,null.min()-.02),max(.25,null.max()+.02));ax.set_xlabel('Spearman correlation under shuffled external gene identities',fontsize=9.7);ax.set_ylabel('Permutation draws',fontsize=10);ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,fontsize=8.7,loc='upper left')
fig.text(.12,.95,'TCGA → CPTAC-2 | DNA-union rank concordance',fontsize=13,fontweight='bold');fig.text(.12,.899,f"Observed rho = {obs:.3f} | Exceedances = {int(test['Null_exceedances'])}/10,000 | Empirical P = {float(test['Empirical_P_greater']):.2e}",fontsize=9.2)
fig.text(.12,.06,'Endpoint: source putative-driver mutation OR amplification OR deep deletion. A small P value does not imply a large correlation.',fontsize=8.1)
finalize(fig,'Figure_6_Permutation_null_CPTAC_DNA_union','The actual 10,000 null correlations for test '+test['Test_ID']+' are shown. The vertical dashed line is the observed Spearman statistic. The empirical upper-tail probability is (b+1)/(10,001); no P=0 is reported. This global test is not a per-gene validation statistic.')

# Rank stability: display top30, resampling all466, not reselecting for plotting.
fig=plt.figure(figsize=(9.3,9.6));ax=fig.add_axes([.15,.13,.63,.77]);rr=primary[:30]
mid=np.array([float(r['DNA_rank_boot_median']) for r in rr]);lo=np.array([float(r['DNA_rank_boot_CI95_low']) for r in rr]);hi=np.array([float(r['DNA_rank_boot_CI95_high']) for r in rr]);y=np.arange(30)
ax.errorbar(mid,y,xerr=np.vstack([mid-lo,hi-mid]),fmt='o',markersize=3.3,linewidth=.8,capsize=2,label='Bootstrap median and 95% interval')
ax.set_yticks(y);ax.set_yticklabels(top,fontstyle='italic',fontsize=8.8);ax.invert_yaxis();ax.set_xlim(0,max(hi)*1.05);ax.set_xlabel('DNA alteration rank among all 466 candidates',fontsize=10);ax.spines[['right','top','left']].set_visible(False);ax.tick_params(axis='y',length=0);ax.grid(axis='x',linewidth=.4,alpha=.25);ax.set_axisbelow(True)
for i,r in enumerate(rr):ax.text(1.04,i,f"{100*float(r['DNA_top30_boot_frequency']):.1f}%",transform=ax.get_yaxis_transform(),fontsize=8.4,va='center')
ax.text(1.04,-.95,'Top-30 inclusion',transform=ax.get_yaxis_transform(),fontsize=8.4,fontweight='bold');ax.legend(loc='lower right',frameon=False,fontsize=8.5)
fig.text(.15,.96,'Primary ranking | Patient-bootstrap stability',fontsize=13,fontweight='bold');fig.text(.15,.927,'2,000 patient resamples, stratified by colon and rectum; all 466 genes are reranked in every draw.',fontsize=9)
fig.text(.15,.055,'Displayed genes are the original TCGA top 30. Equal observed frequencies retain average ranks; display ties are alphabetical.',fontsize=8.2)
finalize(fig,'Figure_7_Primary_rank_bootstrap','Patient-bootstrap rank stability for the frozen primary top-30 display. Every replicate resamples whole patient columns within colon and rectum and reranks all 466 genes. Intervals are percentile descriptions of rank stability, not simultaneous confidence guarantees or clinical validation.')
(OUT/'FIGURE_LEGENDS.txt').write_text('\n\n'.join(legends),encoding='utf-8')
print('ALL 7 SVG FIGURES SAVED',flush=True)
