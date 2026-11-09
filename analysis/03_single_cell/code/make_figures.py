#!/usr/bin/env python3
"""Regenerate Figure 3A-D from released numerical source CSVs."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def save(fig, outdir: Path, stem: str):
    fig.tight_layout()
    fig.savefig(outdir/(stem+'.svg'), bbox_inches='tight')
    fig.savefig(outdir/(stem+'.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig3a(src: Path, out: Path):
    d=pd.read_csv(src/'Figure_3A_source.csv')
    fig,ax=plt.subplots(figsize=(8,7))
    categories=[
        ('Other candidates', ~d['Spatial_followup_category'].isin(['CORE_47','TIER_A_COMPOSITION_SENSITIVE_18']), 'o', 22),
        ('Tier A, composition-sensitive', d['Spatial_followup_category'].eq('TIER_A_COMPOSITION_SENSITIVE_18'), 'o', 38),
        ('Core 47', d['Spatial_followup_category'].eq('CORE_47'), 'o', 44),
    ]
    for label,mask,marker,size in categories:
        ax.scatter(d.loc[mask,'Tumor_vs_normal_epi_median_delta_log2CPM'],
                   d.loc[mask,'Tumor_epi_vs_nonEpi_median_delta_log2CPM'],
                   label=label,marker=marker,s=size,alpha=.75)
    ax.axhline(0,ls='--',lw=1); ax.axvline(0,ls='--',lw=1)
    ax.set_xlabel('Tumor vs adjacent-normal epithelium\nmedian paired Δlog2(CPM+1)')
    ax.set_ylabel('Tumor epithelium vs tumor non-epithelium\nmedian paired Δlog2(CPM+1)')
    ax.set_title('Single-cell prioritization of 466 COREAD dependency candidates')
    ax.legend(frameon=False)
    for gene in ['SOX9','YAP1','HMGA1','MYC','HSPD1','DDX56','BCL2L1','CTNNB1']:
        r=d.loc[d.Gene.eq(gene)]
        if len(r):
            x=float(r.iloc[0]['Tumor_vs_normal_epi_median_delta_log2CPM']); y=float(r.iloc[0]['Tumor_epi_vs_nonEpi_median_delta_log2CPM'])
            ax.annotate(gene,(x,y),xytext=(5,5),textcoords='offset points',fontsize=8)
    save(fig,out,'Figure_3A_singlecell_priority_landscape')


def fig3b(src: Path, out: Path):
    d=pd.read_csv(src/'Figure_3B_source.csv').set_index('Gene')
    # Archived Figure_3B_source.csv already contains log2(mean CPM + 1).
    vals=d.astype(float)
    z=vals.sub(vals.mean(axis=1),axis=0).div(vals.std(axis=1).replace(0,np.nan),axis=0).fillna(0)
    fig,ax=plt.subplots(figsize=(7,11))
    im=ax.imshow(z.values,aspect='auto')
    ax.set_xticks(range(z.shape[1]),z.columns,rotation=45,ha='right')
    ax.set_yticks(range(z.shape[0]),z.index,fontsize=7)
    ax.set_title('Core 47: relative expression across tumor lineages')
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.02); cb.set_label('Row-standardized log2(mean CPM+1)')
    save(fig,out,'Figure_3B_core47_lineage_heatmap')


def fig3c(src: Path, out: Path):
    d=pd.read_csv(src/'Figure_3C_source.csv')
    x=d['MMRd_vs_MMRp_epi_median_delta_log2CPM'].astype(float)
    q=d['MMRd_vs_MMRp_epi_BH_FDR'].astype(float).clip(lower=np.nextafter(0,1))
    y=-np.log10(q)
    sig=q<0.05
    fig,ax=plt.subplots(figsize=(8,6))
    ax.scatter(x[~sig],y[~sig],label='FDR ≥ 0.05',s=22,alpha=.7)
    ax.scatter(x[sig],y[sig],label='FDR < 0.05',s=28,alpha=.8)
    ax.axvline(0,ls='--',lw=1); ax.axhline(-np.log10(.05),ls='--',lw=1)
    ax.set_xlabel('MMRd vs MMRp tumor epithelium\nmedian Δlog2(CPM+1)')
    ax.set_ylabel('−log10(BH FDR)')
    ax.set_title('Mismatch-repair-associated epithelial expression')
    ax.legend(frameon=False)
    top=d.assign(_abs=x.abs(),_y=y).loc[sig].sort_values(['_abs','_y'],ascending=False).head(10)
    for r in top.itertuples(index=False):
        ax.annotate(r.Gene,(r.MMRd_vs_MMRp_epi_median_delta_log2CPM,-np.log10(max(r.MMRd_vs_MMRp_epi_BH_FDR,np.nextafter(0,1)))),xytext=(5,5),textcoords='offset points',fontsize=8)
    save(fig,out,'Figure_3C_MMRd_vs_MMRp_epithelial')


def fig3d(src: Path, out: Path):
    d=pd.read_csv(src/'Figure_3D_source.csv')
    labels={'CORE_47':'Core 47','TIER_A_COMPOSITION_SENSITIVE_18':'Tier A\ncomposition-sensitive','TIER_B_134':'Tier B','TIER_C_267':'Tier C'}
    x=[labels.get(v,v) for v in d.Category]
    y=d.Genes.astype(int).to_numpy()
    fig,ax=plt.subplots(figsize=(8,6))
    bars=ax.bar(x,y)
    for b,n in zip(bars,y): ax.text(b.get_x()+b.get_width()/2,n+4,str(n),ha='center',va='bottom')
    ax.set_ylabel('Genes'); ax.set_title('Evidence-based single-cell prioritization')
    ax.spines[['top','right']].set_visible(False)
    save(fig,out,'Figure_3D_priority_categories')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]); ap.add_argument('--output',type=Path,required=True); args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    src=args.root/'figures'
    fig3a(src,args.output); fig3b(src,args.output); fig3c(src,args.output); fig3d(src,args.output)
    print('Figures written to',args.output)

if __name__=='__main__':main()
