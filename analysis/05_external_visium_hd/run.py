#!/usr/bin/env python3
"""Recompute corrected external Visium HD validation, ST12 and Figure 5.

processed: use verified 466-gene sample-metric checkpoints;
counts: recompute those metrics from included candidate raw counts and COMPLETE
        library-size denominators; raw: first extract counts from public H5.
All modes rebuild set endpoints, exact tumour-level tests, seeded matched nulls
and gene-level source data. The candidate hierarchy is frozen before validation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from utils import csv_write, json_write, bh, flag, compare_frames

ROOT=Path(__file__).resolve().parent
SAMPLES=['P1CRC','P2CRC','P5CRC']
SETCOL={'Core47':'Core47','Spatial43':'Spatially_substantiated_Core47',
        'Convergent9':'Convergent9','Extended13':'Extended13_5pct'}


def summaries(annotations, sample_metrics):
    sets={name:annotations.loc[flag(annotations[col]),'Gene'].tolist() for name,col in SETCOL.items()}
    setrows=[]
    convergent=[]
    territories=[]
    for sample in SAMPLES:
        d=annotations.merge(sample_metrics[sample],on='Gene',validate='one_to_one')
        for name,genes in sets.items():
            sub=d.loc[d.Gene.isin(genes) & flag(d.evaluable)]
            row={'Sample':sample,'Set':name,'N_evaluable':len(sub)}
            for metric,label,prefix in [('tumor_minus_tissue','Tumor_minus_tissue','Tumor'),
                                        ('periphery_minus_interior','Periphery_minus_interior','Periphery'),
                                        ('Moran_I','Moran_I','Moran')]:
                values=sub[metric].dropna().to_numpy(float)
                row[label+'_mean']=np.mean(values)
                row[label+'_median']=np.median(values)
                row[prefix+'_positive_fraction']=np.mean(values>0)
            row['Status']='Recomputed from complete gene-level metrics'
            setrows.append(row)
        for _,r in d.loc[d.Gene.isin(sets['Convergent9'])].iterrows():
            convergent.append({'Sample':sample,'Gene':r.Gene,
                'Discovery_spatial_category':r.Spatial_category,
                'Discovery_tumor_stroma_effect':r.TumorStroma_median_effect,
                'Discovery_Moran_I':r.Moran_median_I,
                'External_tumor_minus_tissue':r.tumor_minus_tissue,
                'External_periphery_minus_interior':r.periphery_minus_interior,
                'External_mean_expression':r.mean_log2CPM1_all,
                'External_detection_fraction':r.detection_fraction_all,
                'External_Moran_I':r.Moran_I,
                'Status':'Complete gene-level external metrics' if r.evaluable else 'Not evaluable on source panel; not a zero or negative result'})
        for category in ['Interface-persistent','Tumor-interior enriched']:
            sub=d.loc[d.Core47 & d.Spatial_category.eq(category) & d.periphery_minus_interior.notna()]
            concordant=sub.periphery_minus_interior.gt(0) if category=='Interface-persistent' else sub.periphery_minus_interior.lt(0)
            territories.append({'Sample':sample,'Category':category,'N_evaluable':len(sub),
                                'N_concordant':int(concordant.sum()),'Concordant_fraction':float(concordant.mean()),
                                'Expected_direction':'positive' if category=='Interface-persistent' else 'negative'})
    setdf=pd.DataFrame(setrows)
    conv=pd.DataFrame(convergent)
    cross=[]
    for gene,g in conv.groupby('Gene',sort=False):
        category=g.Discovery_spatial_category.iloc[0]
        t=g.External_tumor_minus_tissue.dropna()
        m=g.External_Moran_I.dropna()
        p=g.External_periphery_minus_interior.dropna()
        concordant=(p>0) if category=='Interface-persistent' else (p<0)
        n=int(concordant.sum())
        label=('not evaluable' if len(p)<3 else ('replicated' if n==3 else ('partially replicated' if n>=2 else 'not replicated')))
        cross.append({'Gene':gene,'Discovery_spatial_category':category,
            'Tumor_evaluable_samples':len(t),'Tumor_positive_samples':int((t>0).sum()) if len(t) else np.nan,
            'Moran_evaluable_samples':len(m),'Moran_positive_samples':int((m>0).sum()) if len(m) else np.nan,
            'Territory_evaluable_samples':len(p),'Territory_concordant_samples':n if len(p) else np.nan,
            'Territory_expected_direction':'positive' if category=='Interface-persistent' else 'negative',
            'Territory_replication_class':label,
            'Mean_external_tumor_effect':t.mean(),'Median_external_tumor_effect':t.median(),
            'Mean_external_periphery_effect':p.mean(),'Median_external_periphery_effect':p.median(),
            'Mean_external_Moran_I':m.mean(),'Median_external_Moran_I':m.median()})
    tests=[]
    for name in sets:
        g=setdf.loc[setdf.Set.eq(name)]
        for endpoint,alt in [('Tumor_minus_tissue_mean','greater'),('Periphery_minus_interior_mean','less'),('Moran_I_mean','greater')]:
            v=g.sort_values('Sample')[endpoint].to_numpy(float)
            one=stats.wilcoxon(v,alternative=alt,method='exact').pvalue
            two=stats.wilcoxon(v,alternative='two-sided',method='exact').pvalue
            tests.append({'Set':name,'Endpoint':endpoint,'N_samples':len(v),'Mean':np.mean(v),
                          'Median':np.median(v),'Min':min(v),'Max':max(v),'Expected_direction':alt,
                          'Exact_one_sided_P':one,'Exact_two_sided_P':two,'Status':'Complete three-tumour endpoint'})
    testdf=pd.DataFrame(tests)
    testdf['BH_FDR_across_12_one_sided']=bh(testdf.Exact_one_sided_P)
    return setdf,conv,pd.DataFrame(cross),pd.DataFrame(territories),testdf


def pooled_metrics(annotations, sample_metrics):
    columns=['Gene','Core47','Spatially_substantiated_Core47','Convergent9','Extended13_5pct',
             'Spatial_category','Integrated_class','TumorStroma_median_effect','Moran_median_I']
    d=annotations[columns].rename(columns={'Spatially_substantiated_Core47':'Spatial_supported'}).copy()
    for sample,prefix in [('P2CRC','P2'),('P5CRC','P5')]:
        d=d.merge(sample_metrics[sample][['Gene','mean_log2CPM1_all','detection_fraction_all','Moran_I']].rename(
            columns={'mean_log2CPM1_all':prefix+'_mean_expression','detection_fraction_all':prefix+'_detection','Moran_I':prefix+'_Moran_I'}),on='Gene',validate='one_to_one')
    for endpoint in ['mean_expression','detection','Moran_I']:
        d['pooled_'+endpoint]=d[['P2_'+endpoint,'P5_'+endpoint]].mean(axis=1)
    d['both_evaluable']=d[['P2_Moran_I','P5_Moran_I','pooled_mean_expression','pooled_detection']].notna().all(axis=1)
    return d


def matched_nulls(pooled, annotations, n_permutations=10000):
    d=pooled.loc[pooled.both_evaluable].copy().reset_index(drop=True)
    d['expr_bin']=pd.qcut(d.pooled_mean_expression,5,labels=False,duplicates='drop').astype(int)
    d['det_bin']=pd.qcut(d.pooled_detection,5,labels=False,duplicates='drop').astype(int)
    bins=(d.expr_bin*10+d.det_bin).to_numpy(int)
    metric=d.pooled_Moran_I.to_numpy(float)
    coords=d[['pooled_mean_expression','pooled_detection']].to_numpy(float)
    coords=(coords-coords.mean(axis=0))/coords.std(axis=0,ddof=1)
    gene_index={g:i for i,g in enumerate(d.Gene)}
    rows=[]
    nulls={}
    for si,(name,col) in enumerate(SETCOL.items()):
        target=np.array([gene_index[g] for g in annotations.loc[flag(annotations[col]),'Gene'] if g in gene_index],int)
        targets=set(target.tolist())
        background=np.array([i for i in range(len(d)) if i not in targets],int)
        pools={int(b):background[bins[background]==b] for b in np.unique(bins[background])}
        fallback={i:background[np.argsort(((coords[background]-coords[i])**2).sum(axis=1))] for i in target}
        observed=float(metric[target].mean())
        rng=np.random.default_rng(20261011+si)
        null=np.empty(n_permutations,float)
        for k in range(n_permutations):
            available=np.ones(len(d),bool)
            available[target]=False
            chosen=[]
            for i in rng.permutation(target):
                pool=pools.get(int(bins[i]),np.empty(0,int))
                pool=pool[available[pool]]
                if not len(pool):
                    ordered=fallback[i]
                    # Exact archived rule, without replacement within this draw.
                    pool=ordered[available[ordered]][:max(5,min(50,int(available[ordered].sum())))]
                pick=int(pool[rng.integers(len(pool))])
                chosen.append(pick)
                available[pick]=False
            null[k]=metric[np.array(chosen,int)].mean()
        exceed=int(np.sum(null>=observed))
        rows.append({'Set':name,'N_evaluable':len(target),'Observed_mean_pooled_Moran_I':observed,
                     'Null_mean':float(null.mean()),'Null_SD':float(null.std(ddof=1)),
                     'Null_q025':float(np.quantile(null,.025)),'Null_q975':float(np.quantile(null,.975)),
                     'Exceedances_of_10000':exceed,'Empirical_P_greater':(exceed+1)/(n_permutations+1),
                     'Random_seed':20261011+si,'Permutations':n_permutations})
        nulls[name]=null
    return pd.DataFrame(rows),nulls


def concordance(pooled):
    rows=[]
    for comparison,left,right in [('discovery_vs_external','Moran_median_I','pooled_Moran_I'),
                                   ('P2_vs_P5','P2_Moran_I','P5_Moran_I')]:
        for name,col in [('All_evaluable',None),('Core47','Core47'),('Spatial43','Spatial_supported'),('Convergent9','Convergent9')]:
            sub=pooled if col is None else pooled.loc[pooled[col]]
            sub=sub.dropna(subset=[left,right])
            r=stats.spearmanr(sub[left],sub[right])
            rows.append({'Set':name,'Comparison':comparison,'N':len(sub),'Spearman_rho':r.statistic,'P_value':r.pvalue})
    d=pd.DataFrame(rows)
    d['BH_FDR_across_8']=bh(d.P_value)
    return d


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['processed','counts','raw'],default='processed')
    parser.add_argument('--raw-dir',type=Path)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'results')
    parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args()
    out=args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    annotations=pd.read_csv(ROOT/'data/frozen_discovery_annotations.csv')
    checkpoint_dir=ROOT/'data/counts'
    if args.mode=='raw':
        if args.raw_dir is None:parser.error('--raw-dir is required in raw mode')
        from counts import prepare_checkpoint
        checkpoint_dir=out/'count_checkpoints'
        for sample in SAMPLES:
            prepare_checkpoint(sample,args.raw_dir/f'{sample}_filtered_feature_bc_matrix.h5',
                               args.raw_dir/f'{sample}_Metadata.parquet.gz',annotations.Gene.tolist(),checkpoint_dir)
    metrics={}
    qc={}
    for sample in SAMPLES:
        if args.mode in ['counts','raw']:
            from counts import metrics_from_checkpoint
            metrics[sample],qc[sample]=metrics_from_checkpoint(sample,checkpoint_dir)
            csv_write(metrics[sample],out/'sample_metrics'/f'{sample}_metrics.csv')
            json_write(qc[sample],out/'sample_metrics'/f'{sample}_QC.json')
        else:
            metrics[sample]=pd.read_csv(ROOT/'data/metrics'/f'{sample}_metrics.csv')
            qc[sample]=json.loads((ROOT/'data/metrics'/f'{sample}_QC.json').read_text())
        assert len(metrics[sample])==466 and metrics[sample].Gene.nunique()==466
    setdf,conv,cross,territory,tests=summaries(annotations,metrics)
    pooled=pooled_metrics(annotations,metrics)
    matched,nulls=matched_nulls(pooled,annotations)
    concord=concordance(pooled)
    tables={'set_level':setdf,'Convergent9_three_tumours':conv,
            'Convergent9_cross_sample_summary':cross,'territory_concordance':territory,
            'patient_level_tests':tests,'Moran_P2_P5_466':pooled,
            'Moran_matched_permutations_P2_P5':matched,'Moran_concordance':concord}
    keys={'set_level':('Sample','Set'),'Convergent9_three_tumours':('Sample','Gene'),
          'Convergent9_cross_sample_summary':('Gene',),'territory_concordance':('Sample','Category'),
          'patient_level_tests':('Set','Endpoint'),'Moran_P2_P5_466':('Gene',),
          'Moran_matched_permutations_P2_P5':('Set',),'Moran_concordance':('Set','Comparison')}
    comparisons={}
    for name,d in tables.items():
        csv_write(d,out/f'{name}.csv')
        expected=pd.read_csv(ROOT/'reference'/f'{name}.csv')
        comparisons[name]=compare_frames(d,expected,keys[name],rtol=1e-8,atol=1e-10)
    np.savez_compressed(out/'Moran_matched_nulls_P2_P5.npz',**nulls)
    qcrows=[]
    for sample,q in qc.items():
        d=annotations.merge(metrics[sample],on='Gene')
        qcrows.append({'Sample':sample,'Features':q['n_features'],'Bins':q['n_bins'],
                       'Tissue_bins':q['region_counts']['Tissue'],'Tumor_interior_bins':q['region_counts']['Tumor'],
                       'Periphery_50um_bins':q['region_counts']['50 micron'],
                       'Rook_edges':q['rook_edges'],'Evaluable_candidates':int(d.evaluable.sum()),
                       'Evaluable_Core47':int((d.evaluable & d.Core47).sum()),'Source_SHA256':q['source_h5_sha256']})
    qcd=pd.DataFrame(qcrows)
    csv_write(qcd,out/'sample_QC.csv')
    comparisons['sample_QC']=compare_frames(qcd,pd.read_csv(ROOT/'reference/sample_QC.csv')[qcd.columns],('Sample',))
    # Independently archived P2/P5 single-gene spatial metrics are validation targets.
    for sample in ['P2CRC','P5CRC']:
        old=pd.read_csv(ROOT/'reference'/f'legacy_{sample}_moran_466.csv')
        cols=['Gene','mean_log2CPM1_all','detection_fraction_all','Moran_I']
        comparisons['source_vs_legacy_'+sample]=compare_frames(metrics[sample][cols],old[cols],('Gene',),rtol=1e-8,atol=1e-10)
    csv_write(setdf[['Sample','Set','N_evaluable','Tumor_minus_tissue_mean']],out/'source_data/Figure_5A.csv')
    csv_write(setdf[['Sample','Set','N_evaluable','Moran_I_mean']],out/'source_data/Figure_5B.csv')
    csv_write(conv[['Sample','Gene','External_tumor_minus_tissue','Status']],out/'source_data/Figure_5C.csv')
    csv_write(conv[['Sample','Gene','Discovery_spatial_category','External_periphery_minus_interior','Status']],out/'source_data/Figure_5D.csv')
    if not args.no_plots:
        from plot_figure5 import plot_figure5
        plot_figure5(setdf,conv,out/'figures')
    observed=conv.External_tumor_minus_tissue.notna()
    numerical={'total_bins_1595565':int(qcd.Bins.sum())==1595565,
               'P1_corrected_1011093_edges':int(qcd.loc[qcd.Sample.eq('P1CRC'),'Rook_edges'].iloc[0])==1011093,
               '455_evaluable_each_tumour':qcd.Evaluable_candidates.eq(455).all(),
               'C9_24_evaluable_positive_tumour':int(observed.sum())==24 and conv.loc[observed,'External_tumor_minus_tissue'].gt(0).all(),
               'C9_24_positive_Moran':conv.loc[observed,'External_Moran_I'].gt(0).all(),
               'C9_24_negative_boundary':conv.loc[observed,'External_periphery_minus_interior'].lt(0).all(),
               'HMGB1_three_unavailable':conv.loc[conv.Gene.eq('HMGB1'),'External_tumor_minus_tissue'].isna().all(),
               '12_patient_tests_BH_0125':len(tests)==12 and tests.BH_FDR_across_12_one_sided.eq(.125).all()}
    numerical={k:bool(v) for k,v in numerical.items()}
    status={'module':'05_external_visium_hd','mode':args.mode,
            'passed':bool(all(v['passed'] for v in comparisons.values()) and all(numerical.values())),
            'reference_comparisons':comparisons,'numerical_checks':numerical,
            'count_checkpoint_scope':'1,595,565 source bins; 18,085-feature library denominators; 466 candidate count columns',
            'matched_nulls':'40,000 newly recomputed draws; P2/P5 only; fixed seeds 20261011–20261014',
            'limits':['No original sequencing alignment or source tissue/periphery annotation was rerun. Spatial bins are not biological replicates.']}
    json_write(status,out/'status.json')
    print(json.dumps(status,indent=2))
    if not status['passed']:raise SystemExit(1)


if __name__=='__main__':main()
