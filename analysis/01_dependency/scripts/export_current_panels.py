#!/usr/bin/env python3
"""Regenerate current Figure 1 and panel-level source data from recalculation.

The renderer was reconstructed on 2026-10-31. Panel A uses an explicit exact-
intersection layout; panel D preserves the 90 genes and order extracted from the
original Illustrator figure. Numerical inference remains in the original scripts.
"""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
plt.rcParams.update({'svg.fonttype':'none','svg.hashsalt':'coread-figure1-reconstruction','font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})

def main():
    p=argparse.ArgumentParser();p.add_argument('--numerical',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);src=a.output/'source_data';src.mkdir(exist_ok=True)
    go=pd.read_csv(ROOT/'data/metadata/GO_axis_membership_5880.csv');cols=[c for c in go if c.startswith('All ')];axes=[c[4:] for c in cols]
    mem=go[cols].fillna(0).astype(int);membership=pd.concat([go[['GENES']],mem],axis=1);membership.to_csv(src/'Figure1A_membership.csv',index=False)
    pat=mem.dot(2**np.arange(5));counts=pat.value_counts();rows=[]
    for v in range(1,32):
        rows.append({'pattern':v,'intersection':' & '.join(axes[j] for j in range(5) if v&(1<<j)),'genes':int(counts.get(v,0)),**{axes[j]:int(bool(v&(1<<j))) for j in range(5)}})
    intersections=pd.DataFrame(rows).sort_values('genes',ascending=False,kind='stable');intersections.to_csv(src/'Figure1A_exact_intersections.csv',index=False)
    assert len(go)==5880 and (mem.sum(axis=1)>=2).sum()==1720
    cl=pd.read_csv(a.numerical/'cell_lines/cell_line_results_1165.csv');org=pd.read_csv(a.numerical/'organoids/organoid_results_1720.csv');it=pd.read_csv(a.numerical/'cell_lines/integration_results_1165.csv')
    cl=cl[cl.Data_status.eq('EVALUABLE')].copy();org=org[org.Data_status.eq('EVALUABLE')].copy();joint=it[it.Both_evaluable.eq('YES')].copy()
    assert len(cl)==1620 and len(org)==1636 and len(joint)==1613
    cl[['Gene_symbol','N_COREAD_evaluable','N_COREAD_positive','COREAD_dep_frequency','COREAD_P_mean','CL_candidate','Frequency_preference','Fisher_BH_FDR','N_other_evaluable','N_other_positive']].to_csv(src/'Figure1B_source.csv',index=False)
    org[['Gene_symbol','N_3DO_GE','GE_lt_minus0p5_n_3DO','GE_lt_minus0p5_frequency_3DO','Chronos_mean_3DO','Primary_candidate_3DO','Concordant_GE_EB_candidate_3DO','Same_model_GE_EB_n_3DO']].to_csv(src/'Figure1C_source.csv',index=False)
    integration_columns=['Gene_symbol','Evidence_pattern','Union_candidate','Both_systems_candidate','CL_candidate','Organoid_3DO_candidate','Both_evaluable','CL_N','CL_positive','CL_frequency','Org_N_3DO','Org_positive_3DO','Org_frequency_3DO']
    it[integration_columns].to_csv(src/'Figure1E_all_1720.csv',index=False);joint[integration_columns].to_csv(src/'Figure1E_displayed_1613.csv',index=False)
    order=(ROOT/'data/metadata/Figure1D_original_gene_order.txt').read_text().split();assert len(order)==90
    hs=org.set_index('Gene_symbol').loc[order];modelcols=[c for c in hs if c.startswith('Chronos_3DO|')]
    heat=hs[modelcols].rename(columns=lambda c:c.split('|',1)[1]);heat.to_csv(src/'Figure1D_exact_90_profiles.csv')
    hs[['GE_lt_minus0p5_n_3DO','GE_lt_minus0p5_frequency_3DO','Chronos_mean_3DO','Same_model_GE_EB_n_3DO','Concordant_GE_EB_candidate_3DO']].to_csv(src/'Figure1D_exact_90_annotations.csv')
    fig=plt.figure(figsize=(15.5,13.5),layout='constrained');grid=fig.add_gridspec(3,2,width_ratios=[1.1,1],height_ratios=[1,1,1])
    ag=grid[0,0].subgridspec(2,1,height_ratios=[3,1],hspace=.04);ax=fig.add_subplot(ag[0]);dots=fig.add_subplot(ag[1],sharex=ax)
    xx=np.arange(31);ax.bar(xx,intersections.genes,color='#477C9C');ax.set_ylabel('Genes in exact intersection');ax.set_title('A  Multiaxis gene universe',loc='left',fontweight='bold');ax.tick_params(axis='x',bottom=False,labelbottom=False)
    for x,n in zip(xx,intersections.genes):ax.text(x,n+12,str(n),ha='center',va='bottom',fontsize=6)
    for j,axis in enumerate(axes):
        dots.scatter(xx,np.full(31,j),s=9,color='#dce1e6');on=intersections[axis].to_numpy().astype(bool);dots.scatter(xx[on],np.full(on.sum(),j),s=15,color='#243B53')
    dots.set_yticks(range(5),axes,fontsize=7);dots.set_xticks([]);dots.invert_yaxis();dots.spines[['left','bottom']].set_visible(False)
    ax.text(.99,.97,'5,880 genes; 1,720 in ≥2 axes',transform=ax.transAxes,ha='right',va='top',fontsize=8)
    ax=fig.add_subplot(grid[0,1]);cand=cl.CL_candidate.eq('YES');higher=cl.Frequency_preference.eq('HIGHER_IN_COREAD')
    for mask,color,label,marker in [(~cand,'#B6BEC8','Other measured genes','o'),(cand&~higher,'#4C91B5','Dependency candidates','o'),(higher,'#C9514E','Higher frequency in COREAD','D')]:
        ax.scatter(cl.loc[mask,'COREAD_dep_frequency']*100,cl.loc[mask,'COREAD_P_mean'],s=16,alpha=.65,c=color,marker=marker,label=f'{label} (n={mask.sum():,})',linewidths=0)
    ax.axvline(200/63,ls='--',color='#888888',lw=.7);ax.set(xlim=(-3,103),ylim=(-.03,1.06),xlabel='Positive COREAD cell lines (%)',ylabel='Mean dependency probability');ax.set_title('B  Cell-line dependency evidence',loc='left',fontweight='bold');ax.legend(fontsize=7,frameon=False,loc='lower right')
    for gene in ['CTNNB1','KRAS','SOX9','FASN','HNF1A']:
        r=cl.loc[cl.Gene_symbol.eq(gene)].iloc[0];ax.annotate(gene,(100*r.COREAD_dep_frequency,r.COREAD_P_mean),xytext=({'CTNNB1':(4,11),'KRAS':(4,-10),'SOX9':(4,6),'FASN':(4,-8)}.get(gene,(4,4))),textcoords='offset points',fontsize=7)
    ax=fig.add_subplot(grid[1,0]);cand=org.Primary_candidate_3DO.eq('YES');con=org.Concordant_GE_EB_candidate_3DO.eq('YES')
    for mask,color,label in [(~cand,'#B6BEC8','Other measured genes'),(cand&~con,'#D6A34B','Chronos criterion'),(con,'#347D8F','Chronos + fitted mapping concordance')]:
        ax.scatter(org.loc[mask,'GE_lt_minus0p5_frequency_3DO']*100,org.loc[mask,'Chronos_mean_3DO'],s=17,alpha=.7,c=color,label=f'{label} (n={mask.sum():,})',linewidths=0)
    ax.axvline(12.5,ls='--',color='#888888',lw=.7);ax.axhline(-.5,ls='--',color='#888888',lw=.7);ax.set(xlim=(-3,103),xlabel='COREAD 3D organoids with Chronos < −0.5 (%)',ylabel='Mean Chronos gene effect');ax.set_title('C  Organoid dependency evidence',loc='left',fontweight='bold');ax.legend(fontsize=7,frameon=False,loc='lower left')
    ax=fig.add_subplot(grid[1:,1]);im=ax.imshow(heat.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-3,vmax=.5);ax.set_yticks(range(90),heat.index,fontsize=6.5);ax.set_xticks(range(16),heat.columns,rotation=90,fontsize=6.5);ax.set_title('D  Recurrent organoid dependencies (90 genes)',loc='left',fontweight='bold');fig.colorbar(im,ax=ax,fraction=.035,pad=.025,label='Chronos gene effect (display clipped)')
    ax=fig.add_subplot(grid[2,0]);groups=[('NEITHER_BOTH_EVALUABLE','#B6BEC8','Neither criterion','o'),('CL_ONLY_BOTH_EVALUABLE','#6199BF','Cell-line criterion','o'),('ORG_ONLY_BOTH_EVALUABLE','#D6A34B','Organoid criterion','s'),('BOTH_CANDIDATES','#A64055','Both criteria','D')]
    for key,color,label,marker in groups:
        mask=joint.Evidence_pattern.eq(key);ax.scatter(joint.loc[mask,'CL_frequency']*100,joint.loc[mask,'Org_frequency_3DO']*100,c=color,marker=marker,s=19,alpha=.65,label=f'{label} (n={mask.sum():,})',linewidths=0)
    ax.axvline(200/63,color='#888888',lw=.7,ls='--');ax.axhline(12.5,color='#888888',lw=.7,ls='--');ax.set(xlim=(-3,103),ylim=(-3,105),xlabel='Positive COREAD cell lines (%)',ylabel='Positive COREAD 3D organoids (%)');ax.set_title('E  Cross-system dependency recurrence',loc='left',fontweight='bold');ax.legend(fontsize=7,frameon=False,loc='lower right')
    fig.savefig(a.output/'Figure1.svg',metadata={'Date':None});fig.savefig(a.output/'Figure1.png',dpi=180);plt.close(fig)
    (a.output/'figure1_reconstruction.json').write_text(json.dumps({'passed':True,'source':'recomputed numerical outputs','panelA':'exact intersections from original ST_2 membership; presentation reconstructed','panelD':'90 original genes and order verified against Illustrator text','joint_displayed':len(joint),'full_candidate_union':int(it.Union_candidate.eq('YES').sum())},indent=2))

if __name__=='__main__':main()
