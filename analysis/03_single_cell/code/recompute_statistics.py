#!/usr/bin/env python3
"""Reconstruct patient-level inference and compare it with frozen results.

Written 2026-10-31 from the documented analysis specification. Frozen outputs
are used only after calculation, as independent regression references. Original
bootstrap intervals remain archival; the reproducible new intervals are labelled.
"""
from pathlib import Path
import argparse, io, json, time
import numpy as np
import pandas as pd
from scipy import stats
from prepare_pseudobulk import atomic

ROOT=Path(__file__).resolve().parents[1]
LINEAGES=['Epi','TNKILC','Myeloid','B','Plasma','Strom','Mast']

def bh(p):
    p=np.asarray(p,float);q=np.full(len(p),np.nan);ok=np.isfinite(p)
    q[ok]=stats.false_discovery_control(p[ok],method='bh');return q

def paired(d):
    p=[];rb=[]
    for col in d.T:
        nz=col[col!=0]
        p.append(float(stats.wilcoxon(col,alternative='two-sided',method='auto').pvalue) if len(nz) else 1.)
        ranks=stats.rankdata(abs(nz));rb.append(float(np.sum(ranks*np.sign(nz))/np.sum(ranks)) if len(nz) else np.nan)
    return {'n':np.repeat(len(d),d.shape[1]),'median':np.median(d,axis=0),'mean':np.mean(d,axis=0),'p':p,'q':bh(p),'rb':rb,'positive_fraction':np.mean(d>0,axis=0)}

def tau(values):
    mx=values.max(axis=0);return np.where(mx>0,(1-values/np.where(mx>0,mx,1)).sum(axis=0)/(len(values)-1),np.nan)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pseudobulk',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True);start=time.time()
    z=np.load(a.pseudobulk/'pseudobulk_counts.npz');genes=z['genes'].astype(str);n=len(genes);assert n==466
    top=pd.read_csv(a.pseudobulk/'top_groups.csv');fine=pd.read_csv(a.pseudobulk/'fine_groups.csv')
    tc=z['top_counts'];tp=z['top_positive'];fc=z['fine_counts'];fp=z['fine_positive']
    cpm=tc/top.library_size.to_numpy()[:,None]*1e6;log=np.log2(cpm+1)
    fcpm=fc/fine.library_size.to_numpy()[:,None]*1e6;flog=np.log2(fcpm+1)
    mapping=pd.read_csv(ROOT/'input/gene_feature_matches_full466.csv').set_index('Gene').loc[genes]
    result=pd.DataFrame({'Gene':genes,'SC_match_method':mapping.Match_method.to_numpy(),'SC_H5_feature_name':mapping.H5_feature_name.to_numpy(),'SC_H5_feature_id':mapping.H5_feature_id.to_numpy()})
    tum=top.SPECIMEN_TYPE.eq('T');epi=top.clTopLevel.eq('Epi')
    ti=np.flatnonzero(tum&epi);ni=np.flatnonzero(~tum&epi)
    for prefix,ii in [('Tumor',ti),('Normal',ni)]:
        result[prefix+'_epi_patient_prevalence_CPM_ge1']=(cpm[ii]>=1).mean(axis=0)
        result[prefix+'_epi_fraction_positive_cells']=tp[ii].sum(axis=0)/top.iloc[ii].n_cells.sum()
    result['Tumor_epi_patient_detection_any_count']=(tc[ti]>0).mean(axis=0)
    tlookup=dict(zip(top.iloc[ti].PID,ti));nlookup=dict(zip(top.iloc[ni].PID,ni));patients=sorted(set(tlookup)&set(nlookup))
    tid=np.array([tlookup[p] for p in patients]);nid=np.array([nlookup[p] for p in patients]);d=log[tid]-log[nid];s=paired(d)
    pre='Tumor_vs_normal_epi_'
    for suffix,key in [('N_pairs','n'),('median_delta_log2CPM','median'),('mean_delta_log2CPM','mean'),('Wilcoxon_P','p'),('BH_FDR','q'),('rank_biserial','rb')]:result[pre+suffix]=s[key]
    result['Tumor_gt_normal_epi_pair_fraction']=s['positive_fraction']
    rng=np.random.default_rng(20261007);boot_indices=rng.integers(0,len(d),size=(2000,len(d)))
    low=[];high=[]
    for col in d.T:
        ci=np.quantile(np.median(col[boot_indices],axis=1),[.025,.975]);low.append(ci[0]);high.append(ci[1])
    result['Reconstructed_boot_CI95_low']=low;result['Reconstructed_boot_CI95_high']=high
    enough=(top.iloc[tid].n_cells.to_numpy()>=100)&(top.iloc[nid].n_cells.to_numpy()>=100);s100=paired(d[enough])
    result['Sensitivity_ge100cells_N_pairs']=s100['n'];result['Sensitivity_ge100cells_median_delta_log2CPM']=s100['median'];result['Sensitivity_ge100cells_BH_FDR']=s100['q']
    balanced=[];balance_meta=[]
    for patient in patients:
        sub=fine.loc[fine.PID.eq(patient)&fine.n_cells.ge(20)];tt=sub.loc[sub.SPECIMEN_TYPE.eq('T')].set_index('epiFine');nn=sub.loc[sub.SPECIMEN_TYPE.eq('N')].set_index('epiFine');shared=sorted(set(tt.index)&set(nn.index))
        if len(shared)>=2:
            balanced.append(flog[tt.loc[shared].group_id.to_numpy()].mean(axis=0)-flog[nn.loc[shared].group_id.to_numpy()].mean(axis=0))
            balance_meta.append({'PID':patient,'shared_clusters':len(shared),'cluster_names':';'.join(shared)})
    balanced=np.array(balanced);sc=paired(balanced)
    for suffix,key in [('N_pairs','n'),('median_delta_log2CPM','median'),('Wilcoxon_P','p'),('BH_FDR','q')]:result['Composition_balanced_'+suffix]=sc[key]
    non_epi=[]
    for ii in ti:
        mask=top.PID.eq(top.iloc[ii].PID)&tum&~epi
        counts=tc[mask].sum(axis=0);lib=top.loc[mask].library_size.sum();non_epi.append(np.log2(counts/lib*1e6+1))
    en=paired(log[ti]-np.array(non_epi));pre='Tumor_epi_vs_nonEpi_'
    for suffix,key in [('N_pairs','n'),('median_delta_log2CPM','median'),('Wilcoxon_P','p'),('BH_FDR','q'),('rank_biserial','rb')]:result[pre+suffix]=en[key]
    result['Tumor_epi_gt_nonEpi_pair_fraction']=en['positive_fraction']
    mmrd=log[ti[top.iloc[ti].MMRStatus.eq('MMRd').to_numpy()]];mmrp=log[ti[top.iloc[ti].MMRStatus.eq('MMRp').to_numpy()]]
    mw=stats.mannwhitneyu(mmrd,mmrp,alternative='two-sided',method='auto',axis=0)
    result['Tumor_MMRd_N']=len(mmrd);result['Tumor_MMRp_N']=len(mmrp)
    result['MMRd_vs_MMRp_epi_median_delta_log2CPM']=np.median(mmrd,axis=0)-np.median(mmrp,axis=0)
    result['MMRd_vs_MMRp_epi_MannWhitney_P']=mw.pvalue;result['MMRd_vs_MMRp_epi_BH_FDR']=bh(mw.pvalue)
    result['MMRd_vs_MMRp_epi_rank_biserial']=2*mw.statistic/(len(mmrd)*len(mmrp))-1
    lineage_means=[];lineage_long=[]
    for lin in LINEAGES:
        ix=np.flatnonzero(tum&top.clTopLevel.eq(lin));mean=cpm[ix].mean(axis=0);prev=(cpm[ix]>=1).mean(axis=0);fraction=tp[ix].sum(axis=0)/top.iloc[ix].n_cells.sum()
        lineage_means.append(mean);result[lin+'_tumor_mean_CPM']=mean;result[lin+'_tumor_patient_prevalence_CPM_ge1']=prev
        for j,g in enumerate(genes):lineage_long.append({'Gene':g,'Lineage':lin,'N_patients':len(ix),'Mean_CPM':mean[j],'Median_CPM':np.median(cpm[ix,j]),'Patient_prevalence_CPM_ge1':prev[j],'Fraction_positive_cells':fraction[j]})
    lineage_means=np.array(lineage_means);result['Tumor_lineage_dominant']=np.array(LINEAGES)[lineage_means.argmax(axis=0)];result['Tumor_lineage_specificity_tau']=tau(lineage_means)
    fine_means=[];fine_names=[];fine_long=[]
    for cl in sorted(fine.epiFine.unique()):
        ix=np.flatnonzero(fine.epiFine.eq(cl)&fine.SPECIMEN_TYPE.eq('T'))
        if not len(ix):continue
        mean=fc[ix].sum(axis=0)/fine.iloc[ix].library_size.sum()*1e6;prev=(fcpm[ix]>=1).mean(axis=0);fraction=fp[ix].sum(axis=0)/fine.iloc[ix].n_cells.sum()
        fine_means.append(mean);fine_names.append(cl)
        for j,g in enumerate(genes):fine_long.append({'Gene':g,'Tumor_epithelial_subcluster':cl,'N_patients':len(ix),'Pooled_CPM':mean[j],'Patient_prevalence_CPM_ge1':prev[j],'Fraction_positive_cells':fraction[j]})
    fine_means=np.array(fine_means);result['Tumor_epi_subcluster_dominant']=np.array(fine_names)[fine_means.argmax(axis=0)];result['Tumor_epi_subcluster_specificity_tau']=tau(fine_means);result.loc[fine_means.max(axis=0)==0,'Tumor_epi_subcluster_dominant']='Undetected'
    integration=pd.read_csv(ROOT/'input/integration_results_1165.csv').rename(columns={'Gene_symbol':'Gene'}).set_index('Gene').loc[genes]
    external_cols=['Evidence_pattern','Both_systems_candidate','CL_candidate','Organoid_3DO_candidate','CL_frequency','Org_frequency_3DO','CL_P_mean','Org_Chronos_mean','Org_GE_EB_concordant','CL_preference','CL_Fisher_BH_FDR']
    for col in external_cols:result[col]=integration[col].to_numpy()
    dna=pd.read_csv(ROOT/'input/primary_466_results.csv').set_index('Gene').loc[genes]
    for col in ['DNA_display_rank','DNA_driver_CNA_frequency','Driver_mutation_frequency','Amplification_frequency','Deep_deletion_frequency','DNA_top30_boot_frequency']:result[col]=dna[col].to_numpy()
    both=result.Both_systems_candidate.eq('YES').to_numpy();tn=(np.array(s['q'])<.05)&(s['median']>0);ee=(np.array(en['q'])<.05)&(en['median']>0);prev=result.Tumor_epi_patient_prevalence_CPM_ge1.ge(.5).to_numpy();comp=(np.array(sc['q'])<.05)&(sc['median']>0)
    tier=np.where(both&tn&ee&prev,'A',np.where(tn&prev,'B','C'));core=(tier=='A')&comp
    yn=lambda x:np.where(x,'YES','NO')
    result['SC_tumor_enriched_vs_normal']=yn(tn);result['SC_epi_enriched_vs_microenvironment']=yn(ee);result['SC_tumor_epi_prevalent_ge50pct']=yn(prev)
    result['SingleCell_priority_tier']=tier;result['SingleCell_evidence_count']=both.astype(int)+tn.astype(int)+ee.astype(int)+prev.astype(int)
    result['Composition_balanced_positive_FDR05']=yn(comp)
    result['Primary_tumor_enrichment_robust_ge100cells']=yn((np.array(s100['q'])<.05)&(s100['median']>0))
    result['Primary_tumor_enrichment_robust_composition']=yn(comp)
    result['Spatial_followup_core_47']=yn(core);result['Spatial_followup_category']=np.where(core,'CORE_47',np.where(tier=='A','TIER_A_COMPOSITION_SENSITIVE_18',np.where(tier=='B','TIER_B_134','TIER_C_267')))
    result=result.sort_values(['SingleCell_priority_tier','SingleCell_evidence_count','Tumor_vs_normal_epi_BH_FDR','Tumor_vs_normal_epi_median_delta_log2CPM'],ascending=[True,False,True,False],kind='stable').reset_index(drop=True)
    result['SingleCell_priority_rank']=np.arange(1,len(result)+1);result['Spatial_followup_core_rank']=np.nan
    iscore=result.Spatial_followup_core_47.eq('YES');result.loc[iscore,'Spatial_followup_core_rank']=np.arange(1,iscore.sum()+1)
    reference=pd.read_csv(ROOT/'results/singlecell_466_results_with_sensitivities.csv').set_index('Gene');actual=result.set_index('Gene').loc[reference.index]
    failures=[];checked=[];maxdiff=0.
    excluded=set(external_cols+['DNA_display_rank','DNA_driver_CNA_frequency','Driver_mutation_frequency','Amplification_frequency','Deep_deletion_frequency','DNA_top30_boot_frequency'])
    for col in actual:
        if col not in reference or col in excluded:continue
        x=actual[col];y=reference[col]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            ok=np.isclose(x,y,rtol=1e-10,atol=1e-10,equal_nan=True);finite=np.isfinite(x)&np.isfinite(y)
            if finite.any():maxdiff=max(maxdiff,float(np.max(abs(x[finite]-y[finite]))))
        else:ok=x.fillna('').astype(str).eq(y.fillna('').astype(str)).to_numpy()
        checked.append(col)
        if not np.all(ok):failures.append({'column':col,'different':int((~ok).sum()),'examples':[{'gene':g,'actual':str(x.loc[g]),'reference':str(y.loc[g])} for g in x.index[~ok][:3]]})
    intervals={}
    for suffix in ['low','high']:
        intervals[suffix+'_max_abs_difference']=float(np.max(abs(actual['Reconstructed_boot_CI95_'+suffix]-reference['Tumor_vs_normal_epi_boot_CI95_'+suffix])))
    counts={'paired_primary':len(d),'paired_ge100cells':int(enough.sum()),'paired_composition_balanced':len(balanced),'tumor_epi_non_epi_pairs':len(ti),'MMRd':len(mmrd),'MMRp':len(mmrp),'positive_primary':int(tn.sum()),'positive_epi_vs_non_epi':int(ee.sum()),'positive_ge100cells':int(((np.array(s100['q'])<.05)&(s100['median']>0)).sum()),'positive_composition_balanced':int(comp.sum()),'MMR_FDR05':int((result.MMRd_vs_MMRp_epi_BH_FDR<.05).sum()),'tiers':result.SingleCell_priority_tier.value_counts().to_dict(),'core47':int(iscore.sum())}
    report={'passed':not failures,'implementation':'reconstructed 2026-10-31','rows':n,'columns_compared':len(checked),'values_compared':n*len(checked),'columns':checked,'max_absolute_difference':maxdiff,'failures':failures,'counts':counts,'bootstrap':{'new_seed':20261007,'draws':2000,'schedule':'one shared patient-index array; lexicographically sorted paired patient IDs; float64 medians; NumPy linear quantiles','historical_schedule_verified':False,'comparison_informational_only':intervals},'seconds':time.time()-start}
    atomic(a.output/'singlecell_466_recomputed.csv',result.to_csv(index=False).encode())
    atomic(a.output/'tumor_lineage_summary_long.csv',pd.DataFrame(lineage_long).to_csv(index=False).encode())
    atomic(a.output/'tumor_epithelial_subcluster_summary_long.csv',pd.DataFrame(fine_long).to_csv(index=False).encode())
    atomic(a.output/'composition_balanced_patient_clusters.csv',pd.DataFrame(balance_meta).to_csv(index=False).encode())
    atomic(a.output/'CORE_47_genes.txt',('\n'.join(result.loc[iscore,'Gene'])+'\n').encode())
    atomic(a.output/'numerical_regression.json',json.dumps(report,indent=2).encode());print(json.dumps(report,indent=2))
    if failures:raise SystemExit('Reconstruction disagrees with frozen reference; inspect numerical_regression.json')

if __name__=='__main__':main()
