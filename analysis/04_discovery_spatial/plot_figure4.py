"""Plot all four current Figure 4 panels from the regenerated ST11 table."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLORS = {'Convergent-9':'#b74436', 'Extended-only sensitivity':'#b88923',
          'Spatially substantiated Core47':'#277f84',
          'Orthogonally anchored Core47 without spatial support':'#71548c',
          'Other Core47':'#5f78a4', 'Tier A composition-sensitive':'#7d94a5',
          'Tier B':'#969da4', 'Tier C':'#b7bcc2'}


def plot_figure4(data, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,
                         'svg.fonttype':'none','axes.spines.top':False,
                         'axes.spines.right':False,'axes.linewidth':.6,
                         'xtick.major.width':.6,'ytick.major.width':.6,
                         'savefig.facecolor':'white'})
    fig = plt.figure(figsize=(12.6,11.2), layout='constrained')
    grid = fig.add_gridspec(2,2, height_ratios=[1,1.35], width_ratios=[1.2,1])
    aa, ab, ac, ad = [fig.add_subplot(grid[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]

    # A. Bubble AREA is directly proportional to the DNA alteration frequency.
    # Missing coordinates stay in the source table but are not plotted.
    a = data.dropna(subset=['CrossSystem_dependency_geomean','Composition_balanced_median_delta_log2CPM'])
    layers=[('Other candidate',~a.Core47,'o','#b7bcc2'),
            ('Core47',a.Core47 & ~a.Convergent9,'o','#277f84'),
            ('Convergent9',a.Convergent9,'D','#b74436')]
    for label,mask,marker,color in layers:
        t=a.loc[mask]
        aa.scatter(t.CrossSystem_dependency_geomean,t.Composition_balanced_median_delta_log2CPM,
                   s=800*t.DNA_driver_CNA_frequency, c=color, marker=marker,
                   alpha=.8 if label!='Other candidate' else .55,
                   linewidths=.55, edgecolors='white', label=label)
    for _,r in a.loc[a.Figure_A_label].iterrows():
        offset={'HSPD1':(-30,-12),'BCL2L1':(5,4),'MYC':(4,5),'FASN':(4,5),'YAP1':(4,-7)}.get(r.Gene,(3,3))
        aa.annotate(r.Gene,(r.CrossSystem_dependency_geomean,r.Composition_balanced_median_delta_log2CPM),
                    xytext=offset,textcoords='offset points',fontsize=6)
    aa.axvline(data.CrossSystem_dependency_geomean.median(),ls='--',lw=.7,color='#91979d')
    aa.axhline(data.Composition_balanced_median_delta_log2CPM.median(),ls='--',lw=.7,color='#91979d')
    aa.set(xlabel='Cross-system dependency (geometric mean)',
           ylabel='Composition-balanced epithelial effect\nΔ log₂(CPM + 1)',title='A   Integrated candidate landscape')
    aa.legend(handles=[Line2D([],[],marker=m,color='none',markerfacecolor=c,markeredgecolor='white',markersize=6,label=l)
                       for l,_,m,c in layers],frameon=False,loc='upper left',fontsize=7)

    b=data.loc[data.Core47]
    cats=[('Tumor-interior enriched','o','#3a7eaf'),('Interface-persistent','s','#248578'),
          ('Interface-composition-associated','^','#b78a32')]
    for label,marker,color in cats:
        t=b.loc[b.Spatial_category.eq(label)]
        ab.scatter(t.TumorStroma_median_effect,t.Moran_median_I,
                   s=110*t.CrossSystem_dependency_geomean,marker=marker,c=color,
                   edgecolors='white',linewidths=.5,label=label)
    t=b.loc[b.Convergent9]
    ab.scatter(t.TumorStroma_median_effect,t.Moran_median_I,
               s=110*t.CrossSystem_dependency_geomean,marker='D',c='#b74436',edgecolors='white',linewidths=.6)
    for _,r in b.loc[b.Figure_B_label].iterrows():
        offset={'MYC':(-31,-7),'SOX9':(4,9),'FASN':(4,9),'DDX56':(-28,-9),
                'RNF6':(-25,-12),'YAP1':(8,5)}.get(r.Gene,(3,3))
        ab.annotate(r.Gene,(r.TumorStroma_median_effect,r.Moran_median_I),xytext=offset,textcoords='offset points',fontsize=6)
    ab.axvline(b.TumorStroma_median_effect.median(),ls='--',lw=.7,color='#91979d')
    ab.axhline(b.Moran_median_I.median(),ls='--',lw=.7,color='#91979d')
    ab.set(xlabel='Tumour − stroma effect, Δ log₂(CPM + 1)',ylabel="Patient median Moran’s I",
           title='B   Spatial architecture of Core47')
    ab.legend(frameon=False,fontsize=6.6,loc='upper right')

    c=data.loc[data.Figure_C_Top24].sort_values('Figure_C_row_order')
    ccols=[col for col in data if col.startswith('Figure_C_scaled_')]
    axes=['Immune_axis','Metabolic_axis','Endocrine_axis','Neural_axis','Stromal_matrix_axis']
    heat=c[ccols].to_numpy(float)
    im=ac.imshow(heat,aspect='auto',cmap='RdBu_r',vmin=-2.5,vmax=2.5,extent=(-.5,6.5,len(c)-.5,-.5))
    binary=c[axes].astype(float).to_numpy()
    ac.imshow(binary,aspect='auto',cmap='Greys',vmin=0,vmax=1,extent=(6.5,11.5,len(c)-.5,-.5))
    ac.axvline(6.5,color='white',lw=2)
    ac.set_xlim(-.5,11.5)
    ac.set_xticks(range(12),['Dependency','Tumour/normal','Epi/non-epi','Composition','Tumour/stroma','Moran’s I','DNA',
                           'Immune','Metabolic','Endocrine','Neural','Stromal–matrix'],rotation=65,ha='right',fontsize=7)
    ac.set_yticks(range(len(c)),c.Gene,fontsize=7)
    ac.set_title('C   Multimodal evidence in the top 24 visual ranks',loc='left')
    for spine in ac.spines.values():spine.set_visible(False)
    ac.tick_params(length=0)
    fig.colorbar(im,ax=ac,orientation='horizontal',shrink=.55,pad=.03,label='Within-panel z-score (continuous evidence)')

    d=data.loc[data.Figure_D_Top25].sort_values('Figure_D_visual_rank')
    y=np.arange(len(d))
    for i,(_,r) in enumerate(d.iterrows()):
        color=COLORS[r.Integrated_class]
        ad.hlines(i,0,r.Integrated_visual_score,color=color,lw=1.3)
        ad.scatter(r.Integrated_visual_score,i,c=color,s=29,marker='D' if r.Convergent9 else 'o',zorder=3)
    ad.set_yticks(y,d.Gene,fontsize=7)
    ad.invert_yaxis()
    ad.set(xlabel='Descriptive integrated visual score',title='D   Top 25 visual ranks')
    ad.grid(axis='x',color='#e6e8ea',lw=.5)
    handles=[Line2D([],[],color=c,marker='D' if name=='Convergent-9' else 'o',markersize=4,lw=1,label=name)
             for name,c in COLORS.items() if name in set(d.Integrated_class)]
    ad.legend(handles=handles,loc='lower right',frameon=False,fontsize=6)
    fig.savefig(output_dir/'Figure_4_recomputed.svg')
    fig.savefig(output_dir/'Figure_4_recomputed.png',dpi=180)
    plt.close(fig)
