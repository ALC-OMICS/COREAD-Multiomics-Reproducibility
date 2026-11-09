#!/usr/bin/env python3
"""Read-only numerical/data-integrity verification of the shipped single-cell release."""
from __future__ import annotations
import argparse, csv, json, math
from collections import Counter
from pathlib import Path

EXPECTED_CORE = 47

def f(x):
    try:return float(x)
    except:return math.nan

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument("--report",type=Path,default=None)
    args=ap.parse_args(); root=args.root
    with (root/'results/singlecell_466_results_with_sensitivities.csv').open(encoding='utf-8',newline='') as fh:
        rows=list(csv.DictReader(fh))
    checks=[]
    def ck(name,cond,detail=""):
        if not cond: raise AssertionError(f"FAIL {name}: {detail}")
        checks.append(f"PASS\t{name}\t{detail}")
    genes=[r['Gene'] for r in rows]
    ck('frozen_universe',len(rows)==466 and len(set(genes))==466,'466 unique genes')
    match=Counter(r['SC_match_method'] for r in rows)
    ck('gene_matching',match.get('Exact_symbol',0)==455 and sum(match.values())-match.get('Exact_symbol',0)==11,str(dict(match)))
    ck('primary_pairs',set(r['Tumor_vs_normal_epi_N_pairs'] for r in rows)=={'36'},'36 paired patients')
    sig=[r for r in rows if f(r['Tumor_vs_normal_epi_BH_FDR'])<0.05]
    pos=[r for r in sig if f(r['Tumor_vs_normal_epi_median_delta_log2CPM'])>0]
    neg=[r for r in sig if f(r['Tumor_vs_normal_epi_median_delta_log2CPM'])<0]
    zero=[r for r in sig if f(r['Tumor_vs_normal_epi_median_delta_log2CPM'])==0]
    ck('primary_direction_counts',(len(sig),len(pos),len(neg),len(zero))==(334,219,112,3),f"sig={len(sig)} pos={len(pos)} neg={len(neg)} zero={len(zero)}")
    ck('ge100_pairs',set(r['Sensitivity_ge100cells_N_pairs'] for r in rows)=={'35'},'35 paired patients')
    ge100=[r for r in rows if f(r['Sensitivity_ge100cells_BH_FDR'])<0.05 and f(r['Sensitivity_ge100cells_median_delta_log2CPM'])>0]
    ck('ge100_positive',len(ge100)==220,'220 positive FDR<0.05')
    ck('composition_pairs',set(r['Composition_balanced_N_pairs'] for r in rows)=={'33'},'33 paired patients')
    comp=[r for r in rows if r['Composition_balanced_positive_FDR05']=='YES']
    ck('composition_positive',len(comp)==198,'198 positive FDR<0.05')
    e_sig=[r for r in rows if f(r['Tumor_epi_vs_nonEpi_BH_FDR'])<0.05]
    e_pos=[r for r in e_sig if f(r['Tumor_epi_vs_nonEpi_median_delta_log2CPM'])>0]
    e_neg=[r for r in e_sig if f(r['Tumor_epi_vs_nonEpi_median_delta_log2CPM'])<0]
    e_zero=[r for r in e_sig if f(r['Tumor_epi_vs_nonEpi_median_delta_log2CPM'])==0]
    ck('epithelium_vs_non_epi',(len(e_sig),len(e_pos),len(e_neg),len(e_zero))==(394,206,179,9),f"sig={len(e_sig)} pos={len(e_pos)} neg={len(e_neg)} zero={len(e_zero)}")
    ck('tumor_epi_vs_non_pairs',set(r['Tumor_epi_vs_nonEpi_N_pairs'] for r in rows)=={'62'},'62 matched tumor patients')
    lineages=Counter(r['Tumor_lineage_dominant'] for r in rows)
    expected={'Strom':140,'Epi':113,'TNKILC':76,'Myeloid':65,'Mast':41,'B':19,'Plasma':12}
    ck('dominant_lineages',dict(lineages)==expected,str(dict(lineages)))
    mmr=[r for r in rows if f(r['MMRd_vs_MMRp_epi_BH_FDR'])<0.05]
    mmr_pos=[r for r in mmr if f(r['MMRd_vs_MMRp_epi_median_delta_log2CPM'])>0]
    mmr_neg=[r for r in mmr if f(r['MMRd_vs_MMRp_epi_median_delta_log2CPM'])<0]
    mmr_zero=[r for r in mmr if f(r['MMRd_vs_MMRp_epi_median_delta_log2CPM'])==0]
    ck('MMR_group_sizes',set(r['Tumor_MMRd_N'] for r in rows)=={'34'} and set(r['Tumor_MMRp_N'] for r in rows)=={'28'},'34 MMRd / 28 MMRp')
    ck('MMR_counts',(len(mmr),len(mmr_pos),len(mmr_neg),len(mmr_zero))==(139,96,40,3),f"sig={len(mmr)} pos={len(mmr_pos)} neg={len(mmr_neg)} zero={len(mmr_zero)}")
    tiers=Counter(r['SingleCell_priority_tier'] for r in rows)
    ck('priority_tiers',dict(tiers)=={'A':65,'B':134,'C':267},str(dict(tiers)))
    core=[r for r in rows if r['Spatial_followup_core_47']=='YES']
    ck('core_47',len(core)==EXPECTED_CORE,'47 genes')
    ck('core_is_tierA',all(r['SingleCell_priority_tier']=='A' for r in core),'all Core 47 are Tier A')
    ck('core_composition',all(r['Composition_balanced_positive_FDR05']=='YES' for r in core),'all Core 47 pass composition-balanced criterion')
    # Numerical source checkpoints, not rendered article graphics.
    for fn,n in [('Figure3A_display_zscores.csv',47),('Figure3B_source.csv',466),('Figure3C_source.csv',4),('Figure3D_source.csv',466)]:
        source=root.parents[1]/'source_data/figure_3'/fn
        with source.open(encoding='utf-8',newline='') as fh: nr=sum(1 for _ in csv.DictReader(fh))
        ck('numeric_source_'+fn,nr==n,f'{nr} rows')
    out='\n'.join(checks)+'\n'
    print(out,end='')
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True); args.report.write_text(out,encoding='utf-8')

if __name__=='__main__':main()
