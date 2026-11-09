"""Patient alteration reanalysis. No synthetic measurements, no inferred driver labels.
Input: immutable OncoPrint exports parsed by parse_exports.py, archived AC-ICAM
public-data summaries, and a frozen 466-gene dependency set.
"""
from pathlib import Path
import csv, json, re, hashlib, math, warnings
from collections import defaultdict, Counter
import numpy as np
from scipy.stats import rankdata, fisher_exact, binomtest, beta, spearmanr
from scipy.special import ndtri
ROOT=Path(__file__).resolve().parents[1]
GENES=(ROOT/'input/union_candidates_466.txt').read_text().split()
G=len(GENES); assert G==466 and len(set(GENES))==G
ENDPOINTS=['DNA_driver_CNA','Driver_mutation','Amplification','Deep_deletion','CNA_any','Any_exported_mutation','DNA_broad','RNA_high','RNA_low','Protein_high','Protein_low']
DNA_ENDPOINTS=ENDPOINTS[:7]
B=10000; BOOT=2000; SEED=20261005
URLS={'COREAD_GDC':'https://www.cbioportal.org/study/summary?id=coad_tcga_gdc; https://www.cbioportal.org/study/summary?id=read_tcga_gdc','COAD_GDC':'https://www.cbioportal.org/study/summary?id=coad_tcga_gdc','READ_GDC':'https://www.cbioportal.org/study/summary?id=read_tcga_gdc','CPTAC2':'https://www.cbioportal.org/study/summary?id=coad_cptac_2019','PANCAN':'https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018','AC_ICAM':'https://github.com/cBioPortal/datahub/tree/master/public/coad_silu_2022'}

def clean(x):
 if isinstance(x,np.ndarray) and x.ndim==0:return clean(x.item())
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple,np.ndarray)):return [clean(v) for v in x]
 if isinstance(x,(np.integer,)):return int(x)
 if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
 if isinstance(x,(np.bool_,)):return bool(x)
 return x

def writecsv(path,rows):
 assert rows
 keys=list(rows[0]);
 for r in rows: assert list(r)==keys, 'Inconsistent output columns'
 with path.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(clean(rows))

def bh(p):
 p=np.array(p,float);out=np.full(p.shape,np.nan);idx=np.flatnonzero(np.isfinite(p))
 if not len(idx):return out
 o=idx[np.argsort(p[idx])];v=p[o]*len(o)/np.arange(1,len(o)+1);out[o]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return out

def wilson(k,n):
 k=np.asarray(k,float);n=np.asarray(n,float)
 with np.errstate(divide='ignore',invalid='ignore'):
  z=ndtri(.975);p=k/n;d=1+z*z/n;mid=(p+z*z/(2*n))/d;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
  lo=np.maximum(0,mid-half);hi=np.minimum(1,mid+half)
 return np.where(n>0,lo,np.nan),np.where(n>0,hi,np.nan)

def frequency(k,n):return np.divide(k,n,out=np.full(np.shape(k),np.nan),where=np.asarray(n)>0)
def ranks(f):
 out=np.full(len(f),np.nan);ok=np.isfinite(f);out[ok]=rankdata(-np.asarray(f)[ok],method='average');return out

def info(key):
 z=dict(np.load(ROOT/f'results/{key}_calls.npz',allow_pickle=False));assert list(z['genes'])==GENES;return z
DATA={k:info(k) for k in ['COAD_GDC','READ_GDC','COREAD_GDC','CPTAC2','PANCAN']}

def get_matrix(key,ep):
 z=DATA[key];shape=(G,len(z['patients']))
 def av(typ):
  raw=z.get('avail__'+typ,np.zeros(shape,bool))
  return np.broadcast_to(raw,shape).copy() & z['track__'+typ][:,None]
 def ev(e):return z['event__'+e].astype(bool)
 if ep=='DNA_driver_CNA':return ev('Driver_mutation')|ev('Amplification')|ev('Deep_deletion'),av('MUTATIONS')&av('CNA')
 if ep=='DNA_strict_explicit':return ev('Driver_explicit')|ev('Amplification')|ev('Deep_deletion'),av('MUTATIONS')&av('CNA')
 if ep=='DNA_broad':return ev('Any_exported_mutation')|ev('Amplification')|ev('Deep_deletion'),av('MUTATIONS')&av('CNA')
 if ep=='CNA_any':return ev('Amplification')|ev('Deep_deletion'),av('CNA')
 typ='MUTATIONS' if ep in ['Driver_mutation','Any_exported_mutation','Driver_explicit','Splice_driver_code'] else 'CNA' if ep in ['Amplification','Deep_deletion'] else 'MRNA' if ep.startswith('RNA') else 'PROTEIN'
 # No gene-level protein coverage in uploads: do not use profile flags as a denominator.
 if typ=='PROTEIN':return ev(ep),np.zeros(shape,bool)
 return ev(ep),av(typ)

SUM={}
for key in DATA:
 SUM[key]={}
 for ep in ENDPOINTS+['DNA_strict_explicit','Driver_explicit','Splice_driver_code']:
  x,m=get_matrix(key,ep);n=m.sum(1);k=(x&m).sum(1);f=frequency(k,n);lo,hi=wilson(k,n)
  SUM[key][ep]={'n':n,'k':k,'f':f,'lo':lo,'hi':hi,'rank':ranks(f),'observed_calls':x.sum(1)}

# Public AC-ICAM values, no alias remapping. Mutation eligibility conservatively
# restricted to exact symbols represented in the CNA annotation universe.
def readcounts(fn):
 return {r['Gene']:r for r in csv.DictReader((ROOT/'downloads'/fn).open())}
cna=readcounts('acicam_cna_counts.csv');mut=readcounts('acicam_mutation_counts.csv');rna=readcounts('acicam_rna_cytoband.csv')
for a in [cna,mut,rna]:assert list(a)==GENES

def num(r,c):return float(r[c]) if r[c] not in ['','NA'] else np.nan
SUM['AC_ICAM']={}
for ep in ENDPOINTS:
 n=np.zeros(G,int);k=np.zeros(G,int)
 for i,g in enumerate(GENES):
  r=cna[g];s=mut[g];t=rna[g]
  if ep in ['Amplification','Deep_deletion','CNA_any','Any_exported_mutation','DNA_broad'] and np.isfinite(num(r,'N')):
   n[i]=int(num(r,'N'))
   if ep=='Amplification':k[i]=int(num(r,'AMP'))
   elif ep=='Deep_deletion':k[i]=int(num(r,'DEL'))
   elif ep=='CNA_any':k[i]=int(num(r,'AMP')+num(r,'DEL'))
   elif ep=='Any_exported_mutation':k[i]=int(num(s,'AnyCodingSplice_mut'))
   elif ep=='DNA_broad':k[i]=int(num(s,'DNA_anycoding_or_CNA'))
  elif ep in ['RNA_high','RNA_low'] and np.isfinite(num(t,'N_RNA')):
   n[i]=int(num(t,'N_RNA'));k[i]=int(num(t,ep))
 f=frequency(k,n);lo,hi=wilson(k,n);assert np.all(k<=n)
 SUM['AC_ICAM'][ep]={'n':n,'k':k,'f':f,'lo':lo,'hi':hi,'rank':ranks(f),'observed_calls':np.where(n>0,k,np.nan)}

# Deterministic primary display order; average ranks quantify tied frequencies.
primary=SUM['COREAD_GDC']['DNA_driver_CNA'];porder=np.array(sorted(range(G),key=lambda i:(-primary['f'][i],GENES[i])))
top=porder[:30];topset=set(top);displayrank=np.empty(G,int);displayrank[porder]=np.arange(1,G+1)
(ROOT/'results/primary_top30_genes.txt').write_text('\n'.join(GENES[i] for i in top)+'\n')

# Patient-bootstrap ranking stability. Resample complete patient columns,
# separately within colon and rectal strata of the joint DNA-profile set.
x,m=get_matrix('COREAD_GDC','DNA_driver_CNA');joint=np.all(m,axis=0); assert np.all(m==joint[None,:])
ids=np.flatnonzero(joint);coadsize=len(DATA['COAD_GDC']['patients']);strata=[ids[ids<coadsize],ids[ids>=coadsize]]
rng=np.random.default_rng(SEED);xb=x[:,ids].astype(float)
weights=np.zeros((BOOT,len(ids)),np.float64);start=0
for ss in strata:
 n=len(ss);weights[:,start:start+n]=rng.multinomial(n,np.full(n,1/n),size=BOOT);start+=n
bootfreq=weights@xb.T/len(ids)
bootrank=rankdata(-bootfreq,axis=1,method='average');ranklow,rankmed,rankhi=np.quantile(bootrank,[.025,.5,.975],axis=0)
alpha=np.argsort(np.array(GENES));bootsorted=alpha[np.argsort(-bootfreq[:,alpha],axis=1,kind='stable')[:,:30]];inclusion=np.bincount(bootsorted.ravel(),minlength=G)/BOOT
np.savez_compressed(ROOT/'results/primary_rank_bootstrap.npz',genes=GENES,rank=bootrank.astype(np.float32),frequencies=bootfreq.astype(np.float32))

# Primary site heterogeneity: all DNA endpoints in ONE family (7 * 466).
siteP=[];sitekeys=[]
for ep in DNA_ENDPOINTS:
 a,b=SUM['COAD_GDC'][ep],SUM['READ_GDC'][ep]
 for i in range(G):
  p=fisher_exact([[a['k'][i],a['n'][i]-a['k'][i]],[b['k'][i],b['n'][i]-b['k'][i]]]).pvalue if min(a['n'][i],b['n'][i])>0 else np.nan
  siteP.append(p);sitekeys.append((ep,i))
site={k:(p,q) for k,p,q in zip(sitekeys,siteP,bh(siteP))}

# External per-gene rate comparisons. Significant results indicate HETEROGENEITY,
# never "validation"; no significance threshold is used to filter the candidate set.
external=[];forfdr=defaultdict(list)
for cohort in ['AC_ICAM','CPTAC2','PANCAN']:
 for i in porder:
  for ep in ENDPOINTS:
   a=SUM['COREAD_GDC'][ep];c=SUM['COAD_GDC'][ep];b=SUM[cohort][ep]
   pa=qa=pc=qc=np.nan;npaired=apos=bpos=aonly=bonly=agreement=np.nan
   reason=''
   if cohort=='PANCAN':
    role='PAIRED_TCGA_ROBUSTNESS_NOT_EXTERNAL'
    if ep in DNA_ENDPOINTS:
     xa,ma=get_matrix('COREAD_GDC',ep);xb2,mb=get_matrix('PANCAN',ep)
     pmap={p:j for j,p in enumerate(DATA['COREAD_GDC']['patients'])};ja=np.array([pmap[p] for p in DATA['PANCAN']['patients']]);mask=ma[i,ja]&mb[i]
     v=xa[i,ja][mask];w=xb2[i][mask];npaired=len(v);apos=int(v.sum());bpos=int(w.sum());aonly=int((v&~w).sum());bonly=int((w&~v).sum());agreement=float((v==w).mean()) if len(v) else np.nan
     pa=binomtest(aonly,aonly+bonly,.5).pvalue if aonly+bonly else 1.
     reason='Paired exact McNemar; same TCGA patients, source pipelines/queries may differ.'
    else:reason='Not comparable: split RNA query has unknown profile denominator; protein gene-level coverage unavailable.'
   else:
    role='INDEPENDENT_COLON_COHORT'
    if ep in DNA_ENDPOINTS and b['n'][i]>0:
     pa=fisher_exact([[a['k'][i],a['n'][i]-a['k'][i]],[b['k'][i],b['n'][i]-b['k'][i]]]).pvalue
     pc=fisher_exact([[c['k'][i],c['n'][i]-c['k'][i]],[b['k'][i],b['n'][i]-b['k'][i]]]).pvalue
     reason='Two-sided Fisher rate heterogeneity, not a test proving replication or equivalence.'
     if cohort=='AC_ICAM' and ep in ['Any_exported_mutation','DNA_broad']:reason+=' Broad coding/splice proxy; driver status unavailable.'
    elif ep.startswith('RNA'):reason='Source-specific RNA calls; reference/threshold not harmonized with TCGA export. Descriptive only.'
    elif cohort=='AC_ICAM' and ep in ['DNA_driver_CNA','Driver_mutation']:reason='Driver annotation missing from retrieved MAF; NOT inferred from gene or variant class.'
    else:reason='Unavailable assay or no exact-symbol match; not a biological zero.'
   status='EVALUABLE_PROFILE_BASED' if b['n'][i]>0 else 'NOT_EVALUABLE'
   if cohort=='AC_ICAM' and b['n'][i]>0:status='EXACT_SYMBOL_PUBLIC_MATRIX'
   if cohort=='PANCAN' and ep not in DNA_ENDPOINTS:status='NO_VERIFIED_GENE_PROFILE_DENOMINATOR'
   # Basic unpooled Newcombe/Wilson difference interval.
   delta=a['f'][i]-b['f'][i]
   dlow=delta-np.sqrt((a['f'][i]-a['lo'][i])**2+(b['hi'][i]-b['f'][i])**2)
   dhi=delta+np.sqrt((a['hi'][i]-a['f'][i])**2+(b['f'][i]-b['lo'][i])**2)
   row={'Gene':GENES[i],'Primary_DNA_display_rank':int(displayrank[i]),'Cohort':cohort,'Cohort_role':role,'Endpoint':ep,'Status':status,'TCGA_n_profiled':int(a['n'][i]),'TCGA_altered':int(a['k'][i]) if a['n'][i]>0 else None,'TCGA_frequency':a['f'][i],'External_n_evaluable':int(b['n'][i]),'External_altered':int(b['k'][i]) if b['n'][i]>0 else None,'External_frequency':b['f'][i],'External_CI95_low':b['lo'][i],'External_CI95_high':b['hi'][i],'External_endpoint_rank':b['rank'][i],'TCGA_minus_external':delta,'Difference_CI95_low':dlow,'Difference_CI95_high':dhi,'Heterogeneity_P':pa,'Heterogeneity_BH_FDR':qa,'COAD_n_profiled':int(c['n'][i]),'COAD_altered':int(c['k'][i]) if c['n'][i]>0 else None,'COAD_frequency':c['f'][i],'COAD_heterogeneity_P':pc,'COAD_heterogeneity_BH_FDR':qc,'Paired_n':npaired,'Paired_TCGA_altered':apos,'Paired_PanCancer_altered':bpos,'TCGA_only_call_n':aonly,'PanCancer_only_call_n':bonly,'Paired_call_agreement':agreement,'Observed_export_calls_no_coverage_claim':b['observed_calls'][i],'Primary_top30':i in topset,'Interpretation':reason,'Source_URL':URLS[cohort]}
   if ep not in DNA_ENDPOINTS:
    row['TCGA_minus_external']=np.nan;row['Difference_CI95_low']=np.nan;row['Difference_CI95_high']=np.nan
   if cohort=='PANCAN' and ep in DNA_ENDPOINTS:
    # Align the displayed comparison denominators with the paired hypothesis test.
    row['TCGA_n_profiled']=npaired;row['TCGA_altered']=apos;row['TCGA_frequency']=apos/npaired if npaired else np.nan
    row['External_n_evaluable']=npaired;row['External_altered']=bpos;row['External_frequency']=bpos/npaired if npaired else np.nan
    row['External_CI95_low'],row['External_CI95_high']=wilson(bpos,npaired)
    row['TCGA_minus_external']=(apos-bpos)/npaired if npaired else np.nan
    # Do not attach an independent-sample CI to paired patient results.
    row['Difference_CI95_low']=np.nan;row['Difference_CI95_high']=np.nan
   j=len(external);external.append(row)
   if np.isfinite(pa):forfdr[cohort+'_vs_COREAD'].append((j,pa,'Heterogeneity_BH_FDR'))
   if np.isfinite(pc):forfdr[cohort+'_vs_COAD'].append((j,pc,'COAD_heterogeneity_BH_FDR'))
for fam,items in forfdr.items():
 for (j,p,col),q in zip(items,bh([it[1] for it in items])):external[j][col]=q

# Primary workbook rows (all 466 retained, DNA rank is descriptive, not drug ranking).
functional={r['Gene_symbol']:r for r in csv.DictReader((ROOT/'input/integration_results_1165.csv').open())}
primaryrows=[]
for i in porder:
 g=GENES[i];fn=functional[g]
 row={'DNA_display_rank':int(displayrank[i]),'Gene':g,'DNA_average_rank':primary['rank'][i],'DNA_rank_boot_median':rankmed[i],'DNA_rank_boot_CI95_low':ranklow[i],'DNA_rank_boot_CI95_high':rankhi[i],'DNA_top30_boot_frequency':inclusion[i],'Functional_pattern':fn['Evidence_pattern'],'Cell_line_candidate':fn['CL_candidate'],'Organoid_candidate':fn['Organoid_3DO_candidate'],'Cell_line_dep_frequency':float(fn['CL_frequency']) if fn['CL_frequency'] else None,'Organoid_dep_frequency':float(fn['Org_frequency_3DO']) if fn['Org_frequency_3DO'] else None,'Organoid_Chronos_mean':float(fn['Org_Chronos_mean']) if fn['Org_Chronos_mean'] else None}
 for ep in ENDPOINTS:
  a=SUM['COREAD_GDC'][ep];row[ep+'_n_profiled']=int(a['n'][i]);row[ep+'_altered']=int(a['k'][i]) if a['n'][i] else None;row[ep+'_frequency']=a['f'][i];row[ep+'_CI95_low']=a['lo'][i];row[ep+'_CI95_high']=a['hi'][i]
 for key in ['COAD_GDC','READ_GDC']:
  a=SUM[key]['DNA_driver_CNA'];row[key+'_DNA_n']=int(a['n'][i]);row[key+'_DNA_altered']=int(a['k'][i]);row[key+'_DNA_frequency']=a['f'][i]
 a=SUM['COREAD_GDC']['DNA_strict_explicit'];row['DNA_explicit_driver_only_frequency']=a['f'][i]
 row['Splice_rec_patient_calls']=int(SUM['COREAD_GDC']['Splice_driver_code']['k'][i])
 x,m=get_matrix('COREAD_GDC','DNA_driver_CNA');s=(DATA['COREAD_GDC']['sample_counts']==1)&m[i];row['Single_sample_only_n']=int(s.sum());row['Single_sample_only_DNA_altered']=int(x[i,s].sum());row['Single_sample_only_DNA_frequency']=float(x[i,s].mean())
 row['Colon_rectum_DNA_heterogeneity_P'],row['Colon_rectum_DNA_BH_FDR']=site[('DNA_driver_CNA',i)]
 row['Coverage_qualification']='Profile flags only; per-gene assay callability not recoverable from categorical export.'
 row['RNA_qualification']='Exported High/Low calls; z threshold not supplied; not a tumor-normal contrast.'
 row['Protein_status']='NOT_EVALUABLE_NO_PROFILE'
 row['Source_URL']=URLS['COREAD_GDC'];primaryrows.append(row)

# Global concordance. Unrestricted label permutation is conditional on this
# selected gene set and does not adjust for gene length or annotation propensity.
# CNA sensitivity uses ARM-AGGREGATED frequencies; genes are not treated as
# independent units in this secondary test.
cyto={g:rna[g]['Cytoband'] for g in GENES};arms={g:(re.match(r'^(\d+|X|Y)[pq]',cyto[g]).group(0) if re.match(r'^(\d+|X|Y)[pq]',cyto[g]) else None) for g in GENES}
def rho(a,b):
 if len(a)<3 or np.ptp(a)==0 or np.ptp(b)==0:return np.nan
 ra=rankdata(a);rb=rankdata(b);ra-=ra.mean();rb-=rb.mean();return float(ra@rb/np.sqrt((ra@ra)*(rb@rb)))

def permtest(a,b,seed):
 a=np.asarray(a);b=np.asarray(b);obs=rho(a,b)
 if not np.isfinite(obs):return obs,np.full(B,np.nan),0,np.nan
 ra=rankdata(a);rb=rankdata(b);ra-=ra.mean();rb-=rb.mean();norm=np.sqrt((ra@ra)*(rb@rb));rg=np.random.default_rng(seed)
 # Independent random order per draw, exact tie-aware statistic; no patient shuffling.
 null=np.empty(B)
 for start in range(0,B,250):
  n=min(250,B-start);idx=np.argsort(rg.random((n,len(a))),axis=1);null[start:start+n]=np.sum(ra[None,:]*rb[idx],axis=1)/norm
 exceed=int(np.sum(null>=obs-1e-12));return obs,null,exceed,(exceed+1)/(B+1)

stats=[];nulls={};statsfamilies=defaultdict(list)
for target in ['AC_ICAM','CPTAC2','PANCAN']:
 for ref in (['COREAD_GDC','COAD_GDC'] if target!='PANCAN' else ['COREAD_GDC']):
  for ep in DNA_ENDPOINTS:
   a=SUM[ref][ep];b=SUM[target][ep];fa=a['f'].copy();fb=b['f'].copy()
   if target=='PANCAN':
    xa,ma=get_matrix(ref,ep);xb2,mb=get_matrix(target,ep);pm={p:j for j,p in enumerate(DATA[ref]['patients'])};ja=np.array([pm[p] for p in DATA[target]['patients']]);mask=ma[:,ja]&mb;n=mask.sum(1);fa=frequency((xa[:,ja]&mask).sum(1),n);fb=frequency((xb2&mask).sum(1),n)
   ok=np.isfinite(fa)&np.isfinite(fb)
   methods=['GENE_LABEL_PERMUTATION']+(['CHROMOSOME_ARM_AGGREGATED_PERMUTATION'] if ep in ['Amplification','Deep_deletion','CNA_any'] else [])
   for method in methods:
    aa=fa[ok];bb=fb[ok];ng=int(ok.sum());unit='gene';nunits=ng;armcounts=0
    if method.startswith('CHROMOSOME'):
     groups=defaultdict(list)
     for i in np.flatnonzero(ok):
      if arms[GENES[i]]:groups[arms[GENES[i]]].append(i)
     groups={k:v for k,v in groups.items() if len(v)>=3};aa=np.array([np.mean(fa[v]) for v in groups.values()]);bb=np.array([np.mean(fb[v]) for v in groups.values()]);unit='chromosome_arm_with_at_least_3_matched_genes';nunits=len(groups);armcounts=sum(map(len,groups.values()))
    seed=SEED+len(stats)*1009
    if ng<3 or nunits<3:obs=np.nan;null=np.full(B,np.nan);exc=0;p=np.nan;status='NOT_EVALUABLE_NO_MATCHED_ENDPOINT'
    else:
     obs,null,exc,p=permtest(aa,bb,seed);status='EVALUATED_CONDITIONAL_CONCORDANCE' if np.isfinite(obs) else 'NOT_EVALUABLE_CONSTANT_PROFILE'
    sid=f'T{len(stats)+1:02d}'
    if np.isfinite(obs):nulls[sid]=null.copy()
    fam='TCGA_REPRESENTATION_ROBUSTNESS' if target=='PANCAN' else 'EXTERNAL_COREAD_CONCORDANCE' if ref=='COREAD_GDC' else 'COLON_ONLY_EXTERNAL_SENSITIVITY'
    note='Positive alignment under random gene identities; not gene-level validation, not functional dependency, not equivalence.'
    if method.startswith('CHROMOSOME'):note='Arm-level means reduce within-arm gene pseudoreplication; between-arm independence not guaranteed; low-power sensitivity.'
    if target=='PANCAN':note='Same TCGA patients with paired coverage; analytical robustness, NOT external validation.'
    if ep in ['Any_exported_mutation','DNA_broad'] and target=='AC_ICAM':note+=' Uses broad coding/splice proxy, not driver-specific calls.'
    mclo,mchi=wilson(exc,B) if np.isfinite(p) else (np.nan,np.nan)
    r={'Test_ID':sid,'Reference':ref,'Comparison':target,'Role':'ROBUSTNESS_NONINDEPENDENT' if target=='PANCAN' else 'EXTERNAL_COLON_COHORT','Endpoint':ep,'Null_scheme':method,'Permutation_unit':unit,'Matched_genes':ng,'Units_tested':nunits,'Genes_in_arm_sensitivity':armcounts,'Spearman_rho':obs,'Permutations':B if np.isfinite(p) else 0,'Seed':seed,'Null_exceedances':exc if np.isfinite(p) else None,'Empirical_P_greater':p,'BH_FDR':np.nan,'BH_family':fam,'Null_mean':float(np.nanmean(null)) if np.isfinite(obs) else np.nan,'Null_SD':float(np.nanstd(null,ddof=1)) if np.isfinite(obs) else np.nan,'Null_percentile_2p5':float(np.nanquantile(null,.025)) if np.isfinite(obs) else np.nan,'Null_percentile_97p5':float(np.nanquantile(null,.975)) if np.isfinite(obs) else np.nan,'Monte_Carlo_tail_CI95_low':mclo,'Monte_Carlo_tail_CI95_high':mchi,'Status':status,'Interpretation':note,'Source_URL':URLS[target]}
    statfamilies_key=fam
    if np.isfinite(p):statsfamilies[fam].append((len(stats),p))
    stats.append(r)
for fam,items in statsfamilies.items():
 for (j,p),q in zip(items,bh([it[1] for it in items])):stats[j]['BH_FDR']=q

# Output patient matrices are preserved in NPZ, event-level calls in CSV.
writecsv(ROOT/'results/primary_466_results.csv',primaryrows)
writecsv(ROOT/'results/external_466_results_long.csv',external)
writecsv(ROOT/'results/permutation_statistics.csv',stats)
writecsv(ROOT/'results/colon_rectum_heterogeneity.csv',[{'Gene':GENES[i],'Endpoint':ep,'Fisher_P':site[(ep,i)][0],'BH_FDR_3262_tests':site[(ep,i)][1]} for ep in DNA_ENDPOINTS for i in range(G)])
# Persist the 10,000 actual null draws for EACH evaluable test, not only P values.
with (ROOT/'results/permutation_null_draws.csv').open('w',newline='') as f:
 w=csv.writer(f);ks=list(nulls);w.writerow(['Permutation_index']+ks)
 for j in range(B):w.writerow([j+1]+[float(nulls[k][j]) for k in ks])
np.savez_compressed(ROOT/'results/permutation_null_draws.npz',**nulls)
summary={'frozen_genes':G,'cohort_patients':{k:len(v['patients']) for k,v in DATA.items()},'profile_counts':{k:{ep:int(v['n'].max()) for ep,v in s.items()} for k,s in SUM.items()},'AC_ICAM_exact_CNA_genes':int((SUM['AC_ICAM']['Amplification']['n']>0).sum()),'AC_ICAM_exact_RNA_genes':int((SUM['AC_ICAM']['RNA_high']['n']>0).sum()),'primary_joint_DNA_patients':len(ids),'primary_joint_DNA_strata':[len(s) for s in strata],'bootstrap_replicates':BOOT,'permutation_replicates':B,'tests_evaluated':len(nulls),'permutation_family_sizes':{k:len(v) for k,v in statsfamilies.items()},'gene_heterogeneity_family_sizes':{k:len(v) for k,v in forfdr.items()},'primary_top30':[GENES[i] for i in top],'top10_details':[{k:r[k] for k in ['Gene','DNA_driver_CNA_altered','DNA_driver_CNA_n_profiled','DNA_driver_CNA_frequency','Driver_mutation_frequency','Amplification_frequency','Deep_deletion_frequency']} for r in primaryrows[:10]],'limits':['Local export denominators are profile based, not gene-callability verified.','Primary TCGA exports aggregate patients including 3 patients with two samples.','AC-ICAM driver annotation missing: no driver-specific external rate or composite.','RNA references/thresholds not harmonized; source-specific calls descriptive only.','Protein coverage unavailable; no formal protein alteration frequencies.','PanCancer is paired TCGA robustness, not independent validation.','Permutation tests describe global correspondence; per-gene Fisher P tests heterogeneity, not proof of validation.']}
(ROOT/'results/analysis_summary.json').write_text(json.dumps(clean(summary),indent=2))
(ROOT/'results/all_endpoint_summaries.json').write_text(json.dumps(clean(SUM)))
print(json.dumps(clean(summary),indent=2))
print('\nEXTERNAL CORE RESULTS')
for r in stats:
 if r['Null_scheme']=='GENE_LABEL_PERMUTATION' and r['Reference']=='COREAD_GDC':print(r['Test_ID'],r['Comparison'],r['Endpoint'],r['Matched_genes'],r['Spearman_rho'],r['Empirical_P_greater'],r['BH_FDR'])
