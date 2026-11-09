#!/usr/bin/env python3
"""Portable raw-discovery route, adapted from the recovered original analysis.

The keyed Moran correction is applied before multiple testing. Fresh null draws
use sorted gene order; they are new draws and do not claim bitwise equivalence
to the historical archived null vectors. Full source and compact-count reconstruction are validated separately from the
fast processed route; see validation/full_source_H5AD_status.json.
"""
import os, re, glob, json, hashlib, math, warnings, shutil, sys
from collections import defaultdict, deque
import numpy as np, pandas as pd, h5py
from scipy import sparse, stats
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
import argparse
from pathlib import Path
parser=argparse.ArgumentParser(description='Discovery spatial analysis from raw H5AD and source annotations; new canonical-order null sampling.')
parser.add_argument('--h5ad',type=Path)
parser.add_argument('--annotations-dir',type=Path)
parser.add_argument('--cell-density',type=Path)
parser.add_argument('--checkpoint-dir',type=Path,help='Use compact source counts, complete libraries, factors and pathology offline')
parser.add_argument('--write-checkpoint-dir',type=Path,help='Save compact reproducible source inputs while reading the H5AD')
parser.add_argument('--upstream-evidence',required=True,type=Path)
parser.add_argument('--output-dir',required=True,type=Path)
args=parser.parse_args()
if args.checkpoint_dir:
    args.annotations_dir=args.annotations_dir or args.checkpoint_dir/'pathology'
    args.cell_density=args.cell_density or args.checkpoint_dir/'W_cell_density_q05.csv.gz'
elif not all([args.h5ad,args.annotations_dir,args.cell_density]):
    parser.error('Provide --h5ad, --annotations-dir and --cell-density, or --checkpoint-dir')
OUT=str(args.output_dir)
os.makedirs(OUT,exist_ok=True)
FIG=os.path.join(OUT,'diagnostic_figures')
os.makedirs(FIG,exist_ok=True)
SEED=20261008
rng=np.random.default_rng(SEED)

# ---------- Load full source H5AD or compact counts with COMPLETE libraries ----------
if args.checkpoint_dir:
    ck=np.load(args.checkpoint_dir/'spot_metadata.npz',allow_pickle=False)
    spot_ids=ck['spot_ids'];symbols=ck['source_feature_symbols'];samples=ck['samples']
    array_row=ck['array_row'];array_col=ck['array_col'];spatial_xy=ck['spatial_xy']
    X=sparse.load_npz(args.checkpoint_dir/'candidate_counts.npz')
    lib=ck['library_sizes'].copy()
    source_feature_count=int(ck['source_feature_count'])
    source_h5_sha256=str(ck['source_h5_sha256'])
else:
    H5=str(args.h5ad)
    h=hashlib.sha256()
    with open(H5,'rb') as handle:
        while chunk:=handle.read(8*1024*1024):h.update(chunk)
    source_h5_sha256=h.hexdigest()
    if source_h5_sha256!='5c839f4e78b4f334e0c19dd3aef173c332b147b8418e854a0ffb6600e05da79a':
        raise ValueError('Source H5AD differs from the recorded analysis version')
    f=h5py.File(H5,'r')
    spot_ids=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in f['obs/spot_id'][:]])
    symbols=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in f['var/SYMBOL'][:]])
    sample_cats=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in f['obs/__categories/sample'][:]])
    samples=sample_cats[f['obs/sample'][:]]
    array_row=f['obs/array_row'][:].astype(int);array_col=f['obs/array_col'][:].astype(int)
    spatial_xy=f['obsm/spatial'][:].astype(float)
    X=sparse.csr_matrix((f['X/data'][:],f['X/indices'][:],f['X/indptr'][:]),shape=tuple(f['X'].attrs['shape']))
    source_feature_count=len(symbols)
    lib=np.asarray(X.sum(axis=1)).ravel().astype(float)
    f.close()
lib[lib<=0]=1.0

# Candidate list/core from prior single-cell analysis
sc=pd.read_csv(args.upstream_evidence)
# Strip derived spatial checkpoints before computing them from raw measurements.
sc=sc[[c for c in sc if not c.startswith(('Adjusted_interface_','Spatial_mean_','Spatial_detection_','Spatial_module'))]]
core_mask=(sc.SingleCell_priority_tier.eq('A') & sc.Composition_balanced_median_delta_log2CPM.gt(0) & sc.Composition_balanced_BH_FDR.lt(.05))
sc['Spatial_followup_core_47']=np.where(core_mask,'YES','NO')
genes=sc['Gene'].astype(str).tolist()
assert len(genes)==466 and len(set(genes))==466
core=set(sc.loc[sc['Spatial_followup_core_47'].astype(str).str.upper().eq('YES'),'Gene'].astype(str))
tierA=set(sc.loc[sc['SingleCell_priority_tier'].astype(str).eq('A'),'Gene'].astype(str))
assert len(core)==47 and len(tierA)==65
sym_to_idx={g:i for i,g in enumerate(symbols)}
alias={'EPRS1':'EPRS','PHB1':'PHB'}
matched_symbol={g:(g if g in sym_to_idx else alias.get(g)) for g in genes}
missing=[g for g in genes if matched_symbol[g] not in sym_to_idx]
if missing: raise RuntimeError(f'Missing genes in spatial matrix: {missing}')
gidx=np.array([sym_to_idx[matched_symbol[g]] for g in genes],dtype=int)
# Candidate raw counts matrix (spots x 466)
C=X[:,gidx].astype(np.float64).tocsr()
if args.write_checkpoint_dir:
    checkpoint=args.write_checkpoint_dir
    checkpoint.mkdir(parents=True,exist_ok=True)
    sparse.save_npz(checkpoint/'candidate_counts.npz',C)
    np.savez_compressed(checkpoint/'spot_metadata.npz',spot_ids=spot_ids.astype(str),
                        source_feature_symbols=np.array([matched_symbol[g] for g in genes],dtype=str),
                        samples=samples.astype(str),array_row=array_row,array_col=array_col,
                        spatial_xy=spatial_xy,library_sizes=lib,source_feature_count=source_feature_count,
                        source_h5_sha256=source_h5_sha256)
    (checkpoint/'pathology').mkdir(exist_ok=True)
    for p in args.annotations_dir.glob('Pathologist_Annotations_*.csv'):
        shutil.copy2(p,checkpoint/'pathology'/p.name)
    pd.read_csv(args.cell_density).to_csv(checkpoint/'W_cell_density_q05.csv.gz',index=False,compression='gzip')
del X
CPM=C.multiply((1e6/lib)[:,None]).toarray()
Y=np.log2(CPM+1.0)
detect=(C.toarray()>0).astype(np.uint8)

# ---------- Metadata / annotations ----------
obs=pd.DataFrame({'spot_id':spot_ids,'sample':samples,'array_row':array_row,'array_col':array_col,'x':spatial_xy[:,0],'y':spatial_xy[:,1]})
obs['patient']=obs['sample'].str.extract(r'(A\d+)')[0]
obs['rep']=obs['sample'].str.extract(r'Rep(\d+)')[0]
obs['barcode']=obs['spot_id'].str.split('_').str[-1]
# handle sample SN123 A938797 without Count prefix exactly by patient+rep mapping
ann_map={}
ann_source={}
for p in sorted(glob.glob(str(args.annotations_dir/'Pathologist_Annotations_*.csv'))):
    bn=os.path.basename(p)
    pat=re.search(r'(A\d+)',bn).group(1)
    rep=re.search(r'Rep(\d+)',bn).group(1)
    df=pd.read_csv(p)
    acol=[c for c in df.columns if 'Pathologist' in c][0]
    for _,r in df.iterrows():
        bc=str(r['Barcode'])
        ann_map[(pat,rep,bc)]=str(r[acol]).strip().lower() if pd.notna(r[acol]) else 'na'
        ann_source[(pat,rep,bc)]=bn
obs['annotation']=[ann_map.get((p,r,b),'unmatched') for p,r,b in zip(obs.patient,obs.rep,obs.barcode)]
if (obs.annotation=='unmatched').any():
    raise RuntimeError(f"Unmatched annotations: {(obs.annotation=='unmatched').sum()}")

def broad(lbl):
    s=lbl.lower().strip()
    if s in ('exclude','na','unmatched'): return 'exclude'
    if s=='tumor': return 'tumor'
    if s.startswith('tumor&stroma'): return 'mixed_tumor_stroma'
    if 'ic aggregate' in s or 'ic aggr' in s: return 'immune_aggregate'
    if s.startswith('stroma') or s.startswith('connective tissue') or s in ('submucosa','lamina propria','muscularis_ic med to high'):
        return 'stroma'
    if 'epithelium' in s or 'epitehlium' in s or s in ('non neo epithelium','glandular tissue'):
        return 'nonneoplastic_epithelium'
    return 'other'
obs['broad']=obs.annotation.map(broad)

# Sample sizes / pathology summary
path_summary=obs.groupby(['patient','sample','broad']).size().unstack(fill_value=0).reset_index()
path_summary.to_csv(os.path.join(OUT,'sample_pathology_counts.csv'),index=False)

# ---------- Hex graph / interface / distances ----------
# Visium array neighbor offsets: same row +/-2 cols, adjacent rows +/-1 col.
offsets=[(0,-2),(0,2),(-1,-1),(-1,1),(1,-1),(1,1)]
obs['interface_tumor']=False; obs['tumor_distance_to_stroma']=np.nan
obs['tumor_interior']=False
adj_by_sample={}
for s,ix_series in obs.groupby('sample').groups.items():
    ix=np.array(list(ix_series),dtype=int)
    lookup={(int(array_row[i]),int(array_col[i])):i for i in ix}
    rows=[]; cols=[]
    neigh_global={i:[] for i in ix}
    for i in ix:
        r,c=int(array_row[i]),int(array_col[i])
        for dr,dc in offsets:
            j=lookup.get((r+dr,c+dc))
            if j is not None:
                neigh_global[i].append(j)
    adj_by_sample[s]=neigh_global
    stromal_sources=[i for i in ix if obs.at[i,'broad'] in ('stroma','mixed_tumor_stroma')]
    # BFS across full tissue graph from stromal/mixed spots
    dist={i:math.inf for i in ix}; dq=deque()
    for i in stromal_sources: dist[i]=0; dq.append(i)
    while dq:
        u=dq.popleft()
        for v in neigh_global[u]:
            if dist[v]>dist[u]+1:
                dist[v]=dist[u]+1; dq.append(v)
    for i in ix:
        if obs.at[i,'broad']=='tumor':
            d=dist[i] if math.isfinite(dist[i]) else np.nan
            obs.at[i,'tumor_distance_to_stroma']=d
            obs.at[i,'interface_tumor']=(d==1)
            obs.at[i,'tumor_interior']=(d>=2 if np.isfinite(d) else True)

# ---------- Per-section and per-patient effects ----------
# helper calculate mean logCPM effect per sample, then average reps per patient
records=[]
for s,ixs in obs.groupby('sample').groups.items():
    ix=np.array(list(ixs),int); pat=obs.loc[ix[0],'patient']
    idx_int=ix[obs.loc[ix,'interface_tumor'].values]
    idx_core=ix[obs.loc[ix,'tumor_interior'].values]
    idx_tum=ix[obs.loc[ix,'broad'].eq('tumor').values]
    idx_str=ix[obs.loc[ix,'broad'].eq('stroma').values]
    for j,g in enumerate(genes):
        rec={'sample':s,'patient':pat,'Gene':g,
             'n_interface':len(idx_int),'n_interior':len(idx_core),'n_tumor':len(idx_tum),'n_stroma':len(idx_str)}
        rec['interface_effect']=np.nan if min(len(idx_int),len(idx_core))<5 else float(Y[idx_int,j].mean()-Y[idx_core,j].mean())
        rec['tumor_vs_stroma_effect']=np.nan if min(len(idx_tum),len(idx_str))<5 else float(Y[idx_tum,j].mean()-Y[idx_str,j].mean())
        # distance gradient in pure tumor, distances 1..5, at least 20 spots and 3 unique distances
        d=obs.loc[idx_tum,'tumor_distance_to_stroma'].to_numpy(float)
        ok=np.isfinite(d); dd=np.minimum(d[ok],5)
        yy=Y[idx_tum[ok],j]
        if ok.sum()>=20 and len(np.unique(dd))>=3:
            rho=stats.spearmanr(dd,yy).statistic # positive = higher inward, negative = interface-high
        else: rho=np.nan
        rec['distance_rho']=rho
        records.append(rec)
sec_effect=pd.DataFrame(records)
sec_effect.to_csv(os.path.join(OUT,'section_gene_effects.csv.gz'),index=False,compression='gzip')
pat_effect=sec_effect.groupby(['patient','Gene'],as_index=False).agg(
    interface_effect=('interface_effect','mean'), tumor_vs_stroma_effect=('tumor_vs_stroma_effect','mean'),
    distance_rho=('distance_rho','mean'), n_sections=('sample','nunique'))
pat_effect.to_csv(os.path.join(OUT,'patient_gene_effects.csv'),index=False)

# gene-level summaries/tests across patient effects
rows=[]
for g in genes:
    d=pat_effect[pat_effect.Gene==g]
    def summarize(col, alternative='two-sided'):
        x=d[col].dropna().to_numpy(float)
        if len(x)<3: return (len(x),np.nan,np.nan,np.nan,np.nan)
        med=float(np.median(x)); mean=float(np.mean(x)); pos=float(np.mean(x>0))
        try: p=float(stats.wilcoxon(x, alternative=alternative, zero_method='wilcox').pvalue)
        except: p=np.nan
        return len(x),med,mean,pos,p
    ni,mi,meani,posi,pi=summarize('interface_effect')
    nt,mt,meant,post,pt=summarize('tumor_vs_stroma_effect')
    nd,md,meand,posd,pdist=summarize('distance_rho')
    rows.append([g,ni,mi,meani,posi,pi,nt,mt,meant,post,pt,nd,md,meand,posd,pdist])
gsum=pd.DataFrame(rows,columns=['Gene','Interface_N_patients','Interface_median_effect','Interface_mean_effect','Interface_positive_patient_fraction','Interface_Wilcoxon_P',
    'TumorStroma_N_patients','TumorStroma_median_effect','TumorStroma_mean_effect','TumorStroma_positive_patient_fraction','TumorStroma_Wilcoxon_P',
    'Distance_N_patients','Distance_median_rho','Distance_mean_rho','Distance_positive_patient_fraction','Distance_Wilcoxon_P'])
for pcol,qcol in [('Interface_Wilcoxon_P','Interface_BH_FDR'),('TumorStroma_Wilcoxon_P','TumorStroma_BH_FDR'),('Distance_Wilcoxon_P','Distance_BH_FDR')]:
    p=gsum[pcol].to_numpy(float); mask=np.isfinite(p); q=np.full(len(p),np.nan); q[mask]=multipletests(p[mask],method='fdr_bh')[1]; gsum[qcol]=q

# ---------- Composition adjustment (tumor spots only) ----------
# Use q05 Cell2Location factors already in obs H5, 40 variables. PCA >=90% variance + section fixed effects.
w=pd.read_csv(args.cell_density)
assert set(w.spot_id)==set(spot_ids)
w=w.set_index('spot_id').loc[spot_ids]
factor_cols=[c for c in w.columns if c.startswith('q05_spot_factors')]
if len(factor_cols)!=40:raise ValueError('The recorded composition model requires the 40 Belgian-reference factors')
Fraw=np.log1p(w[factor_cols].to_numpy(float))
# tumor spots with interface/interior only
tum_idx=np.where(obs.broad.eq('tumor').to_numpy())[0]
Fs=StandardScaler().fit_transform(Fraw[tum_idx])
pca=PCA(n_components=0.90,svd_solver='full',random_state=SEED)
Fp=pca.fit_transform(Fs)
# section dummies, drop first
D=pd.get_dummies(obs.loc[tum_idx,'sample'],drop_first=True,dtype=float).to_numpy()
Z=np.column_stack([np.ones(len(tum_idx)), D, Fp])
# residualize 466 jointly
B=np.linalg.lstsq(Z,Y[tum_idx],rcond=None)[0]
R=Y[tum_idx]-Z@B
# map global->local
loc={g:i for i,g in enumerate(tum_idx)}
comp_recs=[]
for s,ixs in obs.loc[tum_idx].groupby('sample').groups.items():
    ix=np.array(list(ixs),int); pat=obs.loc[ix[0],'patient']
    ai=ix[obs.loc[ix,'interface_tumor'].values]; bi=ix[obs.loc[ix,'tumor_interior'].values]
    if min(len(ai),len(bi))<5: continue
    al=np.array([loc[i] for i in ai]); bl=np.array([loc[i] for i in bi])
    eff=R[al].mean(axis=0)-R[bl].mean(axis=0)
    comp_recs.extend([{'sample':s,'patient':pat,'Gene':g,'Adjusted_interface_effect':float(eff[j])} for j,g in enumerate(genes)])
comp_sec=pd.DataFrame(comp_recs)
comp_pat=comp_sec.groupby(['patient','Gene'],as_index=False).Adjusted_interface_effect.mean()
comp_sec.to_csv(os.path.join(OUT,'composition_adjusted_section_effects.csv'),index=False)
comp_pat.to_csv(os.path.join(OUT,'composition_adjusted_patient_effects.csv'),index=False)
comp_rows=[]
for g in genes:
    x=comp_pat.loc[comp_pat.Gene==g,'Adjusted_interface_effect'].dropna().to_numpy(float)
    if len(x)>=3:
        try:p=stats.wilcoxon(x).pvalue
        except:p=np.nan
        comp_rows.append([g,len(x),np.median(x),np.mean(x),np.mean(x>0),p])
    else: comp_rows.append([g,len(x),np.nan,np.nan,np.nan,np.nan])
comp_sum=pd.DataFrame(comp_rows,columns=['Gene','Adjusted_interface_N_patients','Adjusted_interface_median_effect','Adjusted_interface_mean_effect','Adjusted_interface_positive_patient_fraction','Adjusted_interface_Wilcoxon_P'])
p=comp_sum.Adjusted_interface_Wilcoxon_P.to_numpy(float); mask=np.isfinite(p); q=np.full(len(p),np.nan); q[mask]=multipletests(p[mask],method='fdr_bh')[1]; comp_sum['Adjusted_interface_BH_FDR']=q

# ---------- Moran I per section (observed; spatial organization recurrence) ----------
moran_recs=[]
for s,ixs in obs.groupby('sample').groups.items():
    ix=np.array(list(ixs),int)
    # exclude pathology 'exclude' only; keep all tissue compartments
    ix=ix[~obs.loc[ix,'broad'].eq('exclude').values]
    if len(ix)<20: continue
    local={g:i for i,g in enumerate(ix)}
    rr=[]; cc=[]
    neigh=adj_by_sample[s]
    for gi in ix:
        a=local[gi]
        for gj in neigh[gi]:
            if gj in local:
                rr.append(a);cc.append(local[gj])
    A=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(ix),len(ix)))
    Wsum=A.sum()
    Zm=Y[ix]-Y[ix].mean(axis=0,keepdims=True)
    den=(Zm*Zm).sum(axis=0)
    num=(Zm*(A@Zm)).sum(axis=0)
    I=(len(ix)/Wsum)*num/np.where(den==0,np.nan,den)
    pat=obs.loc[ix[0],'patient']
    for j,g in enumerate(genes): moran_recs.append({'sample':s,'patient':pat,'Gene':g,'Moran_I':float(I[j]),'n_spots':len(ix)})
moran=pd.DataFrame(moran_recs)
moran.to_csv(os.path.join(OUT,'moran_section_466.csv'),index=False)
moran_pat=moran.groupby(['patient','Gene'],as_index=False).Moran_I.mean()
moran_sum=moran_pat.groupby('Gene').Moran_I.agg(['count','median','mean',lambda x:np.mean(x>0)]).reset_index()
moran_sum.columns=['Gene','Moran_N_patients','Moran_median_I','Moran_mean_I','Moran_positive_patient_fraction']
# patient-level signed-rank vs 0 descriptive
mps=[]
for g in genes:
    x=moran_pat.loc[moran_pat.Gene==g,'Moran_I'].to_numpy(float)
    try:p=stats.wilcoxon(x).pvalue if len(x)>=3 else np.nan
    except:p=np.nan
    mps.append(p)
moran_sum=moran_sum.merge(pd.DataFrame({'Gene':genes,'Moran_Wilcoxon_P':mps}),on='Gene',validate='one_to_one')
p=moran_sum.Moran_Wilcoxon_P.to_numpy(float);mask=np.isfinite(p);q=np.full(len(p),np.nan);q[mask]=multipletests(p[mask],method='fdr_bh')[1];moran_sum['Moran_BH_FDR']=q

# ---------- Expression/detectability metrics for matched null ----------
expr=pd.DataFrame({'Gene':genes,'Spatial_mean_logCPM':Y.mean(axis=0),'Spatial_detection_fraction':detect.mean(axis=0), 'Spatial_mean_CPM':CPM.mean(axis=0)})

# Merge main results with prior evidence
res=sc.merge(gsum,on='Gene',how='left').merge(comp_sum,on='Gene',how='left').merge(moran_sum,on='Gene',how='left').merge(expr,on='Gene',how='left')
res['Spatial_Core47']=res.Gene.isin(core)
res['Spatial_TierA65']=res.Gene.isin(tierA)
# categories based on median raw and adjusted interface effects + tumor/stroma
conditions=[]
for _,r in res.iterrows():
    ie=r.Interface_median_effect; ae=r.Adjusted_interface_median_effect; ts=r.TumorStroma_median_effect
    if np.isfinite(ie) and ie>0 and np.isfinite(ae) and ae>0: cat='Interface-persistent'
    elif np.isfinite(ie) and ie>0: cat='Interface-composition-associated'
    elif np.isfinite(ie) and ie<0 and np.isfinite(ts) and ts>0: cat='Tumor-interior enriched'
    elif np.isfinite(ts) and ts<0: cat='Stroma-associated'
    else: cat='No dominant spatial category'
    conditions.append(cat)
res['Spatial_category']=conditions

# ---------- 10,000 matched gene-set permutations ----------
# Stratify background by expression x detection quantile bins; sample one match per core gene when possible.
res_idx=res.set_index('Gene')
background=[g for g in genes if g not in core]
# quantile ranks across all 466, robust bins 5x5
expr_rank=pd.qcut(res['Spatial_mean_logCPM'].rank(method='first'),5,labels=False)
det_rank=pd.qcut(res['Spatial_detection_fraction'].rank(method='first'),5,labels=False)
res['match_bin']=expr_rank.astype(str)+'_'+det_rank.astype(str)
binmap=res.set_index('Gene')['match_bin'].to_dict()
bg_by_bin=defaultdict(list)
for g in background:bg_by_bin[binmap[g]].append(g)
# fallback candidates closest standardized 2d if bin too small
stdx=(res[['Spatial_mean_logCPM','Spatial_detection_fraction']]-res[['Spatial_mean_logCPM','Spatial_detection_fraction']].mean())/res[['Spatial_mean_logCPM','Spatial_detection_fraction']].std()
coord={g:stdx.iloc[i].to_numpy(float) for i,g in enumerate(res.Gene)}
metrics=['Interface_median_effect','Adjusted_interface_median_effect','Moran_median_I']
obs_metrics={m:float(res.loc[res.Spatial_Core47,m].mean()) for m in metrics}
null={m:np.empty(10000) for m in metrics}
for b in range(10000):
    chosen=[]; available=set(background)
    # random order core genes to reduce allocation bias
    for cg in rng.permutation(sorted(core)):
        cand=[g for g in bg_by_bin[binmap[cg]] if g in available]
        if not cand:
            cand=sorted(available)
            d=np.array([np.sum((coord[g]-coord[cg])**2) for g in cand])
            # nearest 20% pool, min 5
            k=min(len(cand),max(5,int(len(cand)*.2)))
            cand=[cand[i] for i in np.argpartition(d,k-1)[:k]]
        pick=rng.choice(cand); chosen.append(pick);available.remove(pick)
    sub=res_idx.loc[chosen]
    for m in metrics:null[m][b]=np.nanmean(sub[m].to_numpy(float))
perm_rows=[]
for m in metrics:
    ob=obs_metrics[m]; arr=null[m]
    # directional: interface/adjusted/Moran higher is hypothesis
    p=(np.sum(arr>=ob)+1)/(len(arr)+1)
    perm_rows.append({'metric':m,'observed_core47_mean':ob,'null_mean':arr.mean(),'null_sd':arr.std(ddof=1),'null_p2.5':np.quantile(arr,.025),'null_p97.5':np.quantile(arr,.975),'empirical_P_greater':p,'exceedances':int(np.sum(arr>=ob))})
perm=pd.DataFrame(perm_rows)
perm.to_csv(os.path.join(OUT,'canonical_order_matched_permutations_summary.csv'),index=False)
np.savez_compressed(os.path.join(OUT,'canonical_order_matched_permutations_nulls.npz'),**null)

# ---------- Core47 spatial modules via average within-section Spearman correlations ----------
core_list=res.loc[res.Spatial_Core47].sort_values('Spatial_followup_core_rank').Gene.tolist()
core_cols=[genes.index(g) for g in core_list]
cor_mats=[]
for s,ixs in obs.groupby('sample').groups.items():
    ix=np.array(list(ixs),int)
    ix=ix[obs.loc[ix,'broad'].eq('tumor').values]
    if len(ix)<50: continue
    mat=Y[ix][:,core_cols]
    corr=pd.DataFrame(mat).corr(method='spearman').to_numpy()
    corr=np.clip(corr,-.999999,.999999)
    cor_mats.append(np.arctanh(corr))
avg_corr=np.tanh(np.nanmean(cor_mats,axis=0))
np.fill_diagonal(avg_corr,1)
dist=1-avg_corr
Zlink=linkage(squareform(dist,checks=False),method='average')
# choose 3 modules as interpretable exploratory summary
mods=fcluster(Zlink,t=3,criterion='maxclust')
module_df=pd.DataFrame({'Gene':core_list,'Spatial_module':mods})
module_df.to_csv(os.path.join(OUT,'core47_spatial_modules.csv'),index=False)
order=leaves_list(Zlink)

# ---------- Patient x core interface matrix ----------
core_pat=pat_effect[pat_effect.Gene.isin(core_list)].pivot(index='Gene',columns='patient',values='interface_effect').reindex(core_list)
core_pat.to_csv(os.path.join(OUT,'core47_patient_interface_matrix.csv'))

# ---------- Figures ----------
# 4A: representative spatial map: pathology + Core47 score
# score = mean gene-wise z across all spots for core47
coreY=Y[:,core_cols]
mu=coreY.mean(axis=0); sd=coreY.std(axis=0); sd[sd==0]=1
score=((coreY-mu)/sd).mean(axis=1)
# choose tumor sample with most pure tumor spots
sample_tum=obs[obs.broad.eq('tumor')].groupby('sample').size().idxmax()
ix=np.array(list(obs.groupby('sample').groups[sample_tum]),int)
fig,axs=plt.subplots(1,2,figsize=(11,5.2))
class_order=['tumor','mixed_tumor_stroma','stroma','nonneoplastic_epithelium','immune_aggregate','other','exclude']
# use matplotlib defaults by separate scatter calls, not manual colors? user didn't ask colors; okay can defaults
for cls in class_order:
    jj=ix[obs.loc[ix,'broad'].eq(cls).values]
    if len(jj): axs[0].scatter(obs.loc[jj,'x'], -obs.loc[jj,'y'], s=8,label=cls,alpha=.75)
axs[0].set_title(f'Pathology | {sample_tum.replace("Count_","")}');axs[0].set_aspect('equal');axs[0].axis('off');axs[0].legend(fontsize=6,loc='best',frameon=False)
scat=axs[1].scatter(obs.loc[ix,'x'], -obs.loc[ix,'y'], c=score[ix],s=9,cmap='coolwarm',vmin=np.nanpercentile(score[ix],2),vmax=np.nanpercentile(score[ix],98))
axs[1].set_title('Core47 spatial score');axs[1].set_aspect('equal');axs[1].axis('off');fig.colorbar(scat,ax=axs[1],shrink=.7,label='Mean z-score')
fig.tight_layout();fig.savefig(os.path.join(FIG,'Diagnostic_4A_spatial_map.svg'),bbox_inches='tight');fig.savefig(os.path.join(FIG,'Diagnostic_4A_spatial_map.png'),dpi=250,bbox_inches='tight');plt.close(fig)

# 4B: landscape all 466
fig,ax=plt.subplots(figsize=(7.4,6.2))
other=res[~res.Spatial_TierA65];ta=res[res.Spatial_TierA65 & ~res.Spatial_Core47];co=res[res.Spatial_Core47]
ax.scatter(other.Interface_median_effect,other.TumorStroma_median_effect,s=18,alpha=.35,label='Other candidates')
ax.scatter(ta.Interface_median_effect,ta.TumorStroma_median_effect,s=32,alpha=.8,label='Tier A, non-Core47')
ax.scatter(co.Interface_median_effect,co.TumorStroma_median_effect,s=42,alpha=.9,label='Core47')
for _,r in co.sort_values('Spatial_followup_core_rank').head(12).iterrows():
    ax.text(r.Interface_median_effect,r.TumorStroma_median_effect,r.Gene,fontsize=7)
ax.axvline(0,lw=.8);ax.axhline(0,lw=.8);ax.set_xlabel('Interface − tumour interior\nmedian patient effect, log2(CPM+1)');ax.set_ylabel('Tumour − stroma\nmedian patient effect, log2(CPM+1)');ax.legend(frameon=False,fontsize=8);ax.set_title('Spatial compartment landscape of 466 COREAD dependency candidates')
fig.tight_layout();fig.savefig(os.path.join(FIG,'Diagnostic_4B_spatial_landscape.svg'),bbox_inches='tight');fig.savefig(os.path.join(FIG,'Diagnostic_4B_spatial_landscape.png'),dpi=250,bbox_inches='tight');plt.close(fig)

# 4C: null distributions for 3 metrics
fig,axs=plt.subplots(1,3,figsize=(12,3.8))
labels={'Interface_median_effect':'Interface effect','Adjusted_interface_median_effect':'Adjusted interface effect','Moran_median_I':"Moran's I"}
for ax,m in zip(axs,metrics):
    ax.hist(null[m],bins=40,alpha=.75)
    ax.axvline(obs_metrics[m],lw=2)
    p=perm.loc[perm.metric==m,'empirical_P_greater'].iloc[0]
    ax.set_title(labels[m]);ax.set_xlabel('Mean across 47 genes');ax.set_ylabel('Matched sets');ax.text(.98,.95,f'Pemp={p:.4g}',transform=ax.transAxes,ha='right',va='top',fontsize=8)
fig.suptitle('Core47 versus 10,000 expression- and detectability-matched 47-gene sets',y=1.02)
fig.tight_layout();fig.savefig(os.path.join(FIG,'Diagnostic_4C_core47_matched_permutations.svg'),bbox_inches='tight');fig.savefig(os.path.join(FIG,'Diagnostic_4C_core47_matched_permutations.png'),dpi=250,bbox_inches='tight');plt.close(fig)

# 4D: heatmap core47 x patients interface effect, row cluster order
mat=core_pat.to_numpy(float)[order]; labels_gene=np.array(core_list)[order]
fig,ax=plt.subplots(figsize=(7.0,10.0))
im=ax.imshow(mat,aspect='auto',cmap='coolwarm',vmin=-np.nanmax(np.abs(mat)),vmax=np.nanmax(np.abs(mat)))
ax.set_xticks(range(mat.shape[1]));ax.set_xticklabels(core_pat.columns,rotation=45,ha='right',fontsize=8)
ax.set_yticks(range(mat.shape[0]));ax.set_yticklabels(labels_gene,fontsize=6.5)
ax.set_title('Core47 interface effects across patients');fig.colorbar(im,ax=ax,shrink=.6,label='Interface − interior log2(CPM+1)')
fig.tight_layout();fig.savefig(os.path.join(FIG,'Diagnostic_4D_core47_patient_heatmap.svg'),bbox_inches='tight');fig.savefig(os.path.join(FIG,'Diagnostic_4D_core47_patient_heatmap.png'),dpi=250,bbox_inches='tight');plt.close(fig)

# Module heatmap supplementary
fig,ax=plt.subplots(figsize=(9,8))
co=avg_corr[np.ix_(order,order)]
im=ax.imshow(co,aspect='equal',cmap='coolwarm',vmin=-1,vmax=1)
lab=np.array(core_list)[order]
ax.set_xticks(range(47));ax.set_xticklabels(lab,rotation=90,fontsize=5);ax.set_yticks(range(47));ax.set_yticklabels(lab,fontsize=5);ax.set_title('Core47 mean within-tumour spatial co-expression')
fig.colorbar(im,ax=ax,shrink=.7,label='Mean Spearman correlation');fig.tight_layout();fig.savefig(os.path.join(FIG,'Supplementary_Core47_spatial_modules.svg'),bbox_inches='tight');plt.close(fig)

# Save result CSV
res=res.merge(module_df,on='Gene',how='left')
res.to_csv(os.path.join(OUT,'Supplementary_Table_11_spatial_466.csv'),index=False)

# QC/provenance summaries
summary={
 'h5_sha256':source_h5_sha256,
 'n_spots':int(len(obs)),'n_genes_matrix':int(source_feature_count),'n_candidates':466,'n_core47':47,'n_tierA65':65,
 'patients':sorted(obs.patient.unique().tolist()),'n_patients':int(obs.patient.nunique()),'n_sections':int(obs['sample'].nunique()),
 'pathology_counts':obs.broad.value_counts().to_dict(),
 'tumor_patients':sorted(obs.loc[obs.broad.eq('tumor'),'patient'].unique().tolist()),
 'n_tumor_patients':int(obs.loc[obs.broad.eq('tumor'),'patient'].nunique()),
 'n_cell2location_factors':len(factor_cols),'cell2location_pca_components_90pct':int(Fp.shape[1]),'cell2location_pca_variance_explained':float(pca.explained_variance_ratio_.sum()),
 'matched_permutations':10000,'seed':SEED,'representative_sample':sample_tum
}
with open(os.path.join(OUT,'QC_summary.json'),'w') as fh: json.dump(summary,fh,indent=2)
# concise text summary
with open(os.path.join(OUT,'RESULTS_SUMMARY.txt'),'w') as fh:
    fh.write(json.dumps(summary,indent=2)+'\n\n')
    fh.write('Matched permutations:\n'+perm.to_string(index=False)+'\n\n')
    fh.write('Core47 spatial categories:\n'+res[res.Spatial_Core47].Spatial_category.value_counts().to_string()+'\n')
    fh.write('\nTop interface effects Core47:\n'+res[res.Spatial_Core47].sort_values('Interface_median_effect',ascending=False)[['Gene','Interface_median_effect','Adjusted_interface_median_effect','TumorStroma_median_effect','Moran_median_I','Spatial_category']].head(20).to_string(index=False)+'\n')

print('DONE')
print(json.dumps(summary,indent=2))
print('\nPERM')
print(perm.to_string(index=False))
print('\nCORE CATS')
print(res[res.Spatial_Core47].Spatial_category.value_counts())
print('\nTOP CORE')
print(res[res.Spatial_Core47].sort_values('Interface_median_effect',ascending=False)[['Gene','Interface_median_effect','Adjusted_interface_median_effect','TumorStroma_median_effect','Moran_median_I','Spatial_category']].head(15).to_string(index=False))

# Final current hierarchy and independent reference validation.
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run import derive_hierarchy,compare_frames,csv_write,json_write
module_root=Path(__file__).resolve().parents[1]
current=derive_hierarchy(res)
current=current.merge(pd.read_csv(module_root/'config/figure_label_selection.csv'),on='Gene',validate='one_to_one')
reference=pd.read_csv(module_root/'reference/discovery_spatial_expected.csv')
current=current[reference.columns]
comparison=compare_frames(current,reference,rtol=1e-8,atol=1e-10)
csv_write(current,Path(OUT)/'ST11_from_source_counts.csv')
status={'module':'04_discovery_spatial','mode':'compact_counts' if args.checkpoint_dir else 'full_source_h5ad',
        'passed':bool(comparison['passed']),'reference_comparison':comparison,
        'source_h5_sha256':source_h5_sha256,'n_spots':int(len(obs)),
        'composition_factors':len(factor_cols),'PCA_components':int(Fp.shape[1]),
        'PCA_variance_explained':float(pca.explained_variance_ratio_.sum()),
        'fresh_null_draws':10000,'fresh_null_order':'sorted; intentionally different from unrecoverable historical Python-set order',
        'scope':'Source counts, full-library normalization, pathology graph, regional effects, composition refit, Moran, hierarchy and exploratory coexpression modules',
        'limits':['Source cell2location factors and source pathology annotations are reused; cell2location and original sequencing alignment are not refitted.',
                  'Historical matched-null vectors remain separate; this route generates new canonical-order draws.']}
json_write(status,Path(OUT)/'status.json')
print(json.dumps(status,indent=2))
if not status['passed']:raise SystemExit(1)
