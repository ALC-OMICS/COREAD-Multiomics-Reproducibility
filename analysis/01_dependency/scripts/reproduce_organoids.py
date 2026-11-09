"""Recompute organoid selection and statistical context from recovered individual values.

This repository reconstruction (2026-10-03) is not the original unavailable
KDE/EM training script. EB probabilities are reconstructed using archived fitted
coefficients, not refitted. The primary candidate rule uses unmodified Chronos.
"""
from pathlib import Path
import os, json, warnings
import numpy as np
from scipy import stats
from scipy.special import expit
from common import ROOT, load_csv, number, write_csv, bh, compare_rows

def run(output=None):
    out=Path(output or os.environ.get('COREAD_OUTPUT_DIR',ROOT/'build'))/'organoids'
    out.mkdir(parents=True,exist_ok=True)
    expected=load_csv(ROOT/'results/csv/organoid_results_1720.csv')
    headers=list(expected[0])
    genes=[r['Gene_symbol'] for r in expected]
    records=load_csv(ROOT/'data/organoids/query_gene_effect_108_screens.csv')
    screens=list(records[0])[1:];matrix_by_gene={r['Gene_symbol']:r for r in records}
    X=np.array([[number(matrix_by_gene.get(g,{}).get(s)) for s in screens] for g in genes])
    allmeta={r['ScreenID']:r for r in load_csv(ROOT/'data/reference/screen_metadata_valores_recuperados.csv')}
    coremeta=load_csv(ROOT/'data/metadata/organoid_coread_models_and_eb_coefficients.csv')
    coreids=[r['ScreenID'] for r in coremeta]; col={s:i for i,s in enumerate(screens)}
    core=np.array([s in coreids for s in screens]);fmt3=np.array([allmeta[s]['ScreenType']=='3DO' for s in screens])
    excluded_ids={'ACH-003377','ACH-003383','ACH-003722'}
    comparator=np.array([not c and allmeta[s]['ModelID'] not in excluded_ids for c,s in zip(core,screens)])
    assert (core.sum(),(core&fmt3).sum(),comparator.sum(),(comparator&fmt3).sum())==(18,16,87,76)
    assert all(allmeta[s]['PassesQC'] in ('1','True','TRUE') for s in screens)
    assert len(set(allmeta[s]['ModelID'] for s in screens))==108
    gh=[h for h in headers if h.startswith('Chronos_3DO|')];gh2=[h for h in headers if h.startswith('Chronos_2DO|')]
    eh=[h for h in headers if h.startswith('EB_REEST_3DO|')];eh2=[h for h in headers if h.startswith('EB_REEST_2DO|')]
    G=X[:,[col[s] for s in coreids]]
    coeff={r['ScreenID']:r for r in coremeta}
    reconstructed=np.column_stack([expit(number(coeff[s]['EB_intercept'])+number(coeff[s]['EB_slope'])*(G[:,j]-number(coeff[s]['EB_centering_constant']))) for j,s in enumerate(coreids)])
    stored=np.array([[number(r.get(next(h for h in eh+eh2 if h.split('|')[1]==s))) for s in coreids] for r in expected])
    maxerr=float(np.nanmax(np.abs(stored-reconstructed)))
    if maxerr>1e-8:raise AssertionError(f'Probability reconstruction failed: {maxerr}')
    # Numerical reconstruction is checked; archived values are retained to preserve precision.
    P=stored; is3=np.array([r['ScreenType']=='3DO' for r in coremeta]);valid=np.isfinite(G).any(1)
    assert valid.sum()==1636
    gen_scores=load_csv(ROOT/'data/organoids/genome_wide_coread_gene_effect_18435x18.csv')
    genome={r['Gene_symbol']:r for r in gen_scores}
    for i,g in enumerate(genes):
        if g in genome:
            direct=np.array([number(genome[g].get(s)) for s in coreids])
            np.testing.assert_allclose(direct,G[i],rtol=0,atol=0,equal_nan=True)
    controls=load_csv(ROOT/'data/metadata/organoid_eb_training_controls.csv')
    assert len(controls)==648 and not (set(r['Gene_symbol'] for r in controls)&set(genes))
    A=X[valid][:,core&fmt3];B=X[valid][:,comparator&fmt3]
    welch=stats.ttest_ind(A,B,axis=1,equal_var=False)
    m1=A.mean(1);m0=B.mean(1);v1=A.var(1,ddof=1)/A.shape[1];v0=B.var(1,ddof=1)/B.shape[1]
    se=np.sqrt(v1+v0);df=(v1+v0)**2/(v1**2/(A.shape[1]-1)+v0**2/(B.shape[1]-1));delta=m1-m0
    ci=stats.t.ppf(.975,df)*se; q=bh(welch.pvalue)
    np.testing.assert_allclose(q,stats.false_discovery_control(welch.pvalue),rtol=1e-12,atol=1e-14)
    aligned_index={i:j for j,i in enumerate(np.flatnonzero(valid))}
    actual=[]
    for i,r0 in enumerate(expected):
        r=dict(r0)
        for mask,label in [(is3,'3DO'),(np.ones(18,bool),'all')]:
            vals=G[i,mask];pvals=P[i,mask];n=int(np.isfinite(vals).sum());k=int((vals<-.5).sum())
            nkey='N_3DO_GE' if label=='3DO' else 'N_all_GE'
            r[nkey]=n
            if n:
                r['GE_lt_minus0p5_n_'+label]=k;r['GE_lt_minus0p5_frequency_'+label]=k/n
                for stat,fun in [('mean',np.nanmean),('median',np.nanmedian),('SD',lambda v:np.nanstd(v,ddof=1)),('min',np.nanmin),('max',np.nanmax)]:r['Chronos_'+stat+'_'+('3DO' if label=='3DO' else 'all')]=float(fun(vals))
                pn=int(np.isfinite(pvals).sum());pk=int((pvals>.5).sum())
                r['N_3DO_EB' if label=='3DO' else 'N_all_EB']=pn
                for stat,fun in [('mean',np.nanmean),('median',np.nanmedian)]:r['EB_probability_'+stat+'_'+label]=float(fun(pvals))
                r['EB_gt0p5_n_'+label]=pk;r['EB_gt0p5_frequency_'+label]=pk/pn
        if valid[i]:
            k=int((G[i,is3]<-.5).sum());both=int(((G[i,is3]<-.5)&(P[i,is3]>.5)).sum())
            r['Primary_candidate_3DO']='YES' if k>=2 else 'NO';r['Expanded_candidate_18']='YES' if (G[i]<-.5).sum()>=2 else 'NO'
            r['Concordant_GE_EB_candidate_3DO']='YES' if both>=2 else 'NO'
            r['Same_model_GE_EB_n_3DO']=both;r['Same_model_GE_EB_frequency_3DO']=both/16
            r['GE_lt_minus0p5_n_2DO']=int((G[i,~is3]<-.5).sum())
            r['Recurrence_3DO']='NO_SCORE_CALL' if k==0 else 'SINGLE_MODEL' if k==1 else 'CANDIDATE_GE80PCT' if k/16>=.8 else 'CANDIDATE_50_TO_LT80PCT' if k/16>=.5 else 'CANDIDATE_25_TO_LT50PCT' if k/16>=.25 else 'CANDIDATE_LT25PCT'
            j=aligned_index[i]
            r.update(Context_N_other_cancer_3DO=76,Context_Chronos_mean_other_3DO=m0[j],Context_Delta_Chronos_3DO=delta[j],Context_Delta_CI95_lower=delta[j]-ci[j],Context_Delta_CI95_upper=delta[j]+ci[j],Context_Welch_P=welch.pvalue[j],Context_Welch_BH_FDR=q[j],Context_preference_FDR05='YES' if q[j]<.05 and delta[j]<0 else 'NO',Context_tests_BH=len(q))
        actual.append(r)
    checked=compare_rows(actual,expected)
    write_csv(out/'organoid_results_1720.csv',headers,actual)
    # Original secondary comparisons, recomputed from the archived full query matrix.
    mw=stats.mannwhitneyu(A,B,axis=1,alternative='two-sided',method='asymptotic',use_continuity=True)
    mwq=bh(mw.pvalue);cliff=2*mw.statistic/(16*76)-1
    inc=core|comparator;Y=X[valid][:,inc].T
    design=np.column_stack([np.ones(105),core[inc].astype(float),(~fmt3[inc]).astype(float)])
    inv=np.linalg.inv(design.T@design);beta=inv@design.T@Y
    residual=Y-design@beta;leverage=np.sum((design@inv)*design,axis=1)
    bread=inv@design.T; hcse=np.sqrt(np.sum(bread[1,:,None]**2*(residual/(1-leverage[:,None]))**2,axis=0))
    hcp=2*stats.t.sf(np.abs(beta[1]/hcse),102);hcq=bh(hcp);hcwidth=stats.t.ppf(.975,102)*hcse
    differential=[]
    for j,i in enumerate(np.flatnonzero(valid)):
        differential.append({'Gene_symbol':genes[i],'Delta_Chronos_3DO':delta[j],'Delta_CI95_lower':delta[j]-ci[j],'Delta_CI95_upper':delta[j]+ci[j],'Welch_t_3DO':welch.statistic[j],'Welch_df_3DO':df[j],'Primary_Welch_P':welch.pvalue[j],'Primary_Welch_BH_FDR':q[j],'Secondary_MWU_P':mw.pvalue[j],'Secondary_MWU_BH_FDR':mwq[j],'Cliff_delta_3DO':cliff[j],'Adjusted_COREAD_beta':beta[1,j],'HC3_SE':hcse[j],'HC3_CI95_lower':beta[1,j]-hcwidth[j],'HC3_CI95_upper':beta[1,j]+hcwidth[j],'Secondary_HC3_P':hcp[j],'Secondary_HC3_BH_FDR':hcq[j]})
    old=load_csv(ROOT/'results/csv/organoid_differential_context.csv')
    checked_context=compare_rows(differential,old,columns=list(differential[0]),rtol=1e-7,atol=1e-10)
    write_csv(out/'organoid_differential_context.csv',list(differential[0]),differential)
    model_records=[]
    for s in screens:
        rec=dict(allmeta[s]);rec['Analysis_group']='COREAD' if s in coreids else ('EXCLUDED' if allmeta[s]['ModelID'] in excluded_ids else 'COMPARATOR');model_records.append(rec)
    write_csv(out/'organoid_108_model_cohort.csv',list(model_records[0]),model_records)
    audit={'evaluable':int(valid.sum()),'primary_candidates':sum(r['Primary_candidate_3DO']=='YES' for r in actual),'concordant_candidates':sum(r['Concordant_GE_EB_candidate_3DO']=='YES' for r in actual),'expanded_candidates':sum(r['Expanded_candidate_18']=='YES' for r in actual),'gene_table_values_checked':checked,'secondary_values_checked':checked_context,'EB_reconstruction_max_absolute_error':maxerr,'EB_refit_performed':False,'original_EB_training_script_recovered':False,'genome_query_values_match':True,'primary_preference':int(((q<.05)&(delta<0)).sum()),'HC3_preferential_genes':[genes[i] for j,i in enumerate(np.flatnonzero(valid)) if hcq[j]<.05 and beta[1,j]<0],'numerical_reconstruction_not_independent_biological_validation':True}
    (out/'verification.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
    return audit

if __name__=='__main__':run()
