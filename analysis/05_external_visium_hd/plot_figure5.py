"""Render the corrected four-panel external-validation figure from source data."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

SAMPLES=['P1CRC','P2CRC','P5CRC']
SETS=['Core47','Spatial43','Convergent9','Extended13']
COLORS={'Core47':'#5f78a4','Spatial43':'#277f84','Convergent9':'#b74436','Extended13':'#b88923'}
MARKERS={'Core47':'o','Spatial43':'s','Convergent9':'D','Extended13':'^'}


def plot_figure5(set_table, convergent, output_dir):
    output_dir=Path(output_dir)
    output_dir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,
                         'savefig.facecolor':'white'})
    fig,axes=plt.subplots(2,2,figsize=(10.8,9.6),layout='constrained',gridspec_kw={'height_ratios':[1,1.4]})
    for ax,endpoint,title,ylabel in [
        (axes[0,0],'Tumor_minus_tissue_mean','A   External tumour enrichment','Set mean tumour − remaining tissue\nΔ log₂(CPM + 1)'),
        (axes[0,1],'Moran_I_mean','B   External spatial organisation','Set mean Moran’s I')]:
        for name in SETS:
            d=set_table.loc[set_table.Set.eq(name)].set_index('Sample').loc[SAMPLES]
            n=int(d.N_evaluable.iloc[0])
            ax.plot(range(3),d[endpoint],color=COLORS[name],marker=MARKERS[name],ms=5,lw=1.1,
                    label=f'{name} (n = {n})')
        ax.set_xticks(range(3),SAMPLES)
        ax.set_xlim(-.17,2.17)
        ax.axhline(0,color='#9ba2a8',lw=.6)
        ax.set_ylabel(ylabel)
        ax.set_title(title,loc='left',fontweight='bold')
        ax.grid(axis='y',color='#e5e8ea',lw=.5)
    axes[0,0].legend(frameon=False,fontsize=7,loc='upper left')
    genes=convergent.loc[convergent.Sample.eq('P1CRC'),'Gene'].tolist()
    for ax,endpoint,title in [
        (axes[1,0],'External_tumor_minus_tissue','C   Convergent9: tumour enrichment'),
        (axes[1,1],'External_periphery_minus_interior','D   Convergent9: boundary effects')]:
        matrix=convergent.pivot(index='Gene',columns='Sample',values=endpoint).reindex(index=genes,columns=SAMPLES)
        values=matrix.to_numpy(float)
        limit=float(np.nanmax(np.abs(values)))
        cmap=plt.get_cmap('RdBu_r').copy()
        cmap.set_bad('#e9ecef')
        im=ax.imshow(values,aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit))
        labels=[]
        for gene in genes:
            if endpoint.endswith('interior'):
                category=convergent.loc[convergent.Gene.eq(gene),'Discovery_spatial_category'].iloc[0]
                tag='interface-persistent' if category=='Interface-persistent' else 'interior-enriched'
                labels.append(f'{gene} [{tag}]')
            else:labels.append(gene)
        ax.set_yticks(range(len(genes)),labels,fontsize=7)
        ax.set_xticks(range(3),SAMPLES)
        ax.tick_params(length=0)
        for i in range(len(genes)):
            for j in range(3):
                value=values[i,j]
                ax.text(j,i,'NA' if not np.isfinite(value) else f'{value:.2f}',ha='center',va='center',
                        fontsize=8,color='white' if np.isfinite(value) and abs(value)>.62*limit else '#232b32')
        ax.set_title(title,loc='left',fontweight='bold',fontsize=9)
        for spine in ax.spines.values():spine.set_visible(False)
        fig.colorbar(im,ax=ax,orientation='horizontal',shrink=.85,pad=.04,
                     label='Δ log₂(CPM + 1)',ticks=[-limit,0,limit],format='%.1f')
    fig.savefig(output_dir/'Figure_5_recomputed.svg')
    fig.savefig(output_dir/'Figure_5_recomputed.png',dpi=180)
    plt.close(fig)
