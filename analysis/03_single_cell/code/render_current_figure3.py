#!/usr/bin/env python3
"""Current Figure 3: one log2 transform of mean patient-level lineage CPM.

Panel order: A lineage heatmap, B priority landscape, C evidence categories,
D MMR contrast. The legacy rendering helper double-transformed an already
logged source table. This renderer begins from recalculated linear mean CPM.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'svg.fonttype':'none','svg.hashsalt':'coread-figure3-corrected','font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
LINEAGES=['Epi','TNKILC','Myeloid','B','Plasma','Strom','Mast']

def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);source=a.output/'source_data';source.mkdir(exist_ok=True)
    d=pd.read_csv(a.results/'singlecell_466_recomputed.csv');core=d.loc[d.Spatial_followup_core_47.eq('YES')].copy();assert len(core)==47 and len(d)==466
    means=core.set_index('Gene')[[x+'_tumor_mean_CPM' for x in LINEAGES]];means.columns=LINEAGES
    logged=np.log2(means+1);z=logged.sub(logged.mean(axis=1),axis=0).div(logged.std(axis=1,ddof=1).replace(0,np.nan),axis=0).fillna(0)
    means.to_csv(source/'Figure3A_mean_patient_CPM.csv');logged.to_csv(source/'Figure3A_log2_mean_CPM.csv');z.to_csv(source/'Figure3A_display_zscores.csv')
    lineage_counts=pd.read_csv(a.results/'tumor_lineage_summary_long.csv')[['Lineage','N_patients']].drop_duplicates();lineage_counts.to_csv(source/'Figure3A_lineage_patient_counts.csv',index=False)
    d.to_csv(source/'Figure3B_source.csv',index=False)
    counts=pd.DataFrame({'Category':['Core47','Other Tier A','Tier B','Tier C'],'Genes':[47,18,134,267]});counts.to_csv(source/'Figure3C_source.csv',index=False)
    mmrcols=['Gene','Tumor_MMRd_N','Tumor_MMRp_N','MMRd_vs_MMRp_epi_median_delta_log2CPM','MMRd_vs_MMRp_epi_MannWhitney_P','MMRd_vs_MMRp_epi_BH_FDR','MMRd_vs_MMRp_epi_rank_biserial','Spatial_followup_core_47'];d[mmrcols].to_csv(source/'Figure3D_source.csv',index=False)
    fig=plt.figure(figsize=(13,12),layout='constrained');grid=fig.add_gridspec(3,2,width_ratios=[.95,1.3],height_ratios=[1.1,.6,1.1])
    ax=fig.add_subplot(grid[:,0]);im=ax.imshow(z,aspect='auto',vmin=-2.3,vmax=2.3,cmap='RdBu_r');ax.set_yticks(range(47),z.index,fontsize=8);ax.set_xticks(range(7),['Epithelial','T/NK/ILC','Myeloid','B cells','Plasma','Stromal','Mast'],rotation=45,ha='right',fontsize=8);ax.set_title('A  Core47 cellular expression',loc='left',fontweight='bold',pad=12);cb=fig.colorbar(im,ax=ax,orientation='horizontal',pad=.035,fraction=.028);cb.set_label('Row z-score of log2(mean patient-level CPM + 1)',fontsize=8)
    ax=fig.add_subplot(grid[0,1]);groups=[('Tier C',d.SingleCell_priority_tier.eq('C'),'#B9C0C9',16),('Tier B',d.SingleCell_priority_tier.eq('B'),'#5596B8',20),('Other Tier A',d.SingleCell_priority_tier.eq('A')&d.Spatial_followup_core_47.ne('YES'),'#D9AD54',27),('Core47',d.Spatial_followup_core_47.eq('YES'),'#A43E58',33)]
    for label,mask,color,size in groups:ax.scatter(d.loc[mask,'Tumor_vs_normal_epi_median_delta_log2CPM'],d.loc[mask,'Tumor_epi_vs_nonEpi_median_delta_log2CPM'],s=size,c=color,alpha=.8,linewidths=0,label=f'{label} (n={mask.sum()})')
    ax.axhline(0,color='#777777',ls='--',lw=.7);ax.axvline(0,color='#777777',ls='--',lw=.7);ax.set_xlabel('Tumour − adjacent-normal epithelium\nMedian paired Δlog2(CPM + 1)');ax.set_ylabel('Tumour epithelial − non-epithelial\nMedian paired Δlog2(CPM + 1)');ax.set_title('B  Patient-level expression evidence',loc='left',fontweight='bold');ax.legend(frameon=False,fontsize=8,loc='lower right')
    for gene in ['HMGA1','MYC','SOX9','YAP1','BCL2L1','CTNNB1']:
        r=d.loc[d.Gene.eq(gene)].iloc[0];ax.annotate(gene,(r.Tumor_vs_normal_epi_median_delta_log2CPM,r.Tumor_epi_vs_nonEpi_median_delta_log2CPM),xytext=(4,4),textcoords='offset points',fontsize=7)
    ax=fig.add_subplot(grid[1,1]);bars=ax.bar(counts.Category,counts.Genes,color=['#A43E58','#D9AD54','#5596B8','#B9C0C9'],width=.65)
    for bar,n in zip(bars,counts.Genes):ax.text(bar.get_x()+bar.get_width()/2,n+4,str(n),ha='center',va='bottom',fontsize=9)
    ax.set_ylim(0,305);ax.set_ylabel('Genes');ax.set_title('C  Fixed single-cell priority categories',loc='left',fontweight='bold')
    ax=fig.add_subplot(grid[2,1]);x=d.MMRd_vs_MMRp_epi_median_delta_log2CPM;q=d.MMRd_vs_MMRp_epi_BH_FDR;y=-np.log10(q.clip(lower=np.finfo(float).tiny));sig=q<.05
    for mask,color,label in [(~sig,'#B9C0C9','BH FDR ≥0.05'),(sig&(x<0),'#407FAD','Higher in MMRp'),(sig&(x>=0),'#A43E58','Higher in MMRd')]:ax.scatter(x[mask],y[mask],c=color,s=22,alpha=.8,linewidths=0,label=f'{label} (n={mask.sum()})')
    ax.axvline(0,color='#777777',ls='--',lw=.7);ax.axhline(-np.log10(.05),color='#777777',ls='--',lw=.7);ax.set_xlabel('MMRd − MMRp tumour epithelium\nDifference of median log2(CPM + 1)');ax.set_ylabel('−log10(BH FDR)');ax.set_title('D  MMR-associated expression (34 MMRd; 28 MMRp)',loc='left',fontweight='bold');ax.legend(frameon=False,fontsize=8,loc='upper left')
    for i in y.nlargest(5).index:
        gene=d.loc[i,'Gene'];offset={'RAC3':(-2,8),'ATP5F1A':(-30,-13),'C1QBP':(4,4),'MAP2K1':(4,4),'ICAM1':(4,4)}.get(gene,(4,4))
        ax.annotate(gene,(x[i],y[i]),xytext=offset,textcoords='offset points',fontsize=7,ha='right' if gene=='RAC3' else 'left')
    ax.margins(y=.1)
    fig.savefig(a.output/'Figure3.svg',metadata={'Date':None});fig.savefig(a.output/'Figure3.png',dpi=180);plt.close(fig)
    (a.output/'figure3_reconstruction.json').write_text(json.dumps({'passed':True,'genes':466,'core':47,'panelA_transform':'log2(mean patient-level CPM + 1), once, followed by row z-score with ddof=1','panel_order':['A lineage','B priority','C categories','D MMR'],'legacy_double_log_corrected':True},indent=2))

if __name__=='__main__':main()
