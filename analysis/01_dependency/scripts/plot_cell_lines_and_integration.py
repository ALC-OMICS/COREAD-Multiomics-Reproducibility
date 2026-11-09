"""Compact repo figures from frozen results. No label jitter, no implicit zero imputation.

These minimal recreations replace the overlapping conversation previews only for
presentation. The archived figures are retained unchanged in archive/.
"""
from pathlib import Path
import os,csv,json
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT,load_csv,number
mpl.rcParams['svg.fonttype']='none'
mpl.rcParams['font.family']='DejaVu Sans'
mpl.rcParams['font.size']=8
mpl.rcParams['axes.linewidth']=.65
mpl.rcParams['svg.hashsalt']='coread-repository-v1'

def save(fig,path,title):
    fig.savefig(path.with_suffix('.svg'),metadata={'Title':title,'Date':None},bbox_inches=None)
    fig.savefig(path.with_suffix('.png'),dpi=180)
    plt.close(fig)

def main():
    out=Path(os.environ.get('COREAD_FIGURE_DIR',ROOT/'build/figures'));out.mkdir(parents=True,exist_ok=True)
    cl=load_csv(ROOT/'results/csv/cell_line_results_1165.csv')
    cl=[r for r in cl if r['Data_status']=='EVALUABLE'];assert len(cl)==1620
    candidates=np.array([r['CL_candidate']=='YES' for r in cl])
    higher=np.array([r['Frequency_preference']=='HIGHER_IN_COREAD' for r in cl])
    x=np.array([number(r['COREAD_dep_frequency'])*100 for r in cl]);y=np.array([number(r['COREAD_P_mean']) for r in cl])
    fig=plt.figure(figsize=(183/25.4,151/25.4));ax=fig.add_axes([.13,.21,.82,.53])
    for mask,marker,alpha,label in [(~candidates,'o',.18,'Not selected (1,285)'),(candidates&~higher,'o',.75,'Other COREAD candidates (312)'),(higher,'D',.9,'Higher frequency in COREAD, BH q < 0.05 (23)')]:
        ax.scatter(x[mask],y[mask],s=12,marker=marker,alpha=alpha,linewidths=0,label=label,rasterized=False)
    ax.set(xlim=(-3,103),ylim=(-.035,1.04),xticks=[0,25,50,75,100],yticks=[0,.25,.5,.75,1])
    ax.set_xlabel('COREAD cell lines with dependency probability > 0.5 (%)',labelpad=7)
    ax.set_ylabel('Mean dependency probability in COREAD',labelpad=7)
    ax.grid(axis='y',alpha=.15,linewidth=.4);ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=7)
    fig.text(.13,.96,'Dependency evidence in COREAD cell lines',fontsize=11,fontweight='bold',va='top')
    fig.text(.13,.91,'1,620 evaluable genes | 63 COREAD models | 335 candidates',fontsize=8,va='top')
    fig.legend(*ax.get_legend_handles_labels(),loc='upper left',bbox_to_anchor=(.12,.874),fontsize=7,frameon=False,markerscale=1.2)
    fig.text(.13,.105,'Candidates: probability > 0.5 in ≥2 COREAD models. Frequencies exclude missing values.',fontsize=6.7)
    fig.text(.13,.077,'Comparative BH q values: Fisher test against 1,099 source-labelled cancer models.',fontsize=6.7)
    fig.text(.13,.049,'Comparative FDR is not an inclusion criterion or the FDR of the candidate catalogue.',fontsize=6.7)
    save(fig,out/'cell_line_dependency_landscape','COREAD cell-line dependency landscape')
    it=load_csv(ROOT/'results/csv/integration_results_1165.csv');p=[r for r in it if r['Both_evaluable']=='YES'];assert len(p)==1613
    x=np.array([number(r['CL_frequency'])*100 for r in p]);y=np.array([number(r['Org_frequency_3DO'])*100 for r in p])
    pattern=np.array([r['Evidence_pattern'] for r in p]);fig=plt.figure(figsize=(183/25.4,165/25.4));ax=fig.add_axes([.13,.23,.82,.50])
    groups=[('NEITHER_BOTH_EVALUABLE','o',.16,'Neither candidate (1,151)'),('CL_ONLY_BOTH_EVALUABLE','o',.75,'Cell-line criterion only (139)'),('ORG_ONLY_BOTH_EVALUABLE','s',.75,'Organoid criterion only (130)'),('BOTH_CANDIDATES','D',.85,'Both criteria (193)')]
    for key,marker,alpha,label in groups:
        mask=pattern==key;ax.scatter(x[mask],y[mask],s=15,marker=marker,alpha=alpha,linewidths=0,label=label,rasterized=False)
    ax.set(xlim=(-3,103),ylim=(-3,103),xticks=[0,25,50,75,100],yticks=[0,25,50,75,100])
    ax.set_xlabel('COREAD cell-line positive-call frequency (%)',labelpad=7)
    ax.set_ylabel('COREAD 3DO positive-call frequency (%)',labelpad=7)
    ax.grid(alpha=.12,linewidth=.4);ax.set_axisbelow(True);ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=7)
    fig.text(.13,.96,'Cell-line and organoid dependency integration',fontsize=11,fontweight='bold',va='top')
    fig.text(.13,.91,'1,613 jointly evaluable genes | 193 shared candidates',fontsize=8,va='top')
    fig.legend(*ax.get_legend_handles_labels(),loc='upper left',bbox_to_anchor=(.12,.876),fontsize=7,frameon=False,markerscale=1.2,ncol=2,columnspacing=1.4)
    fig.text(.13,.143,'Cell lines: probability > 0.5 in ≥2 models; organoids: Chronos < −0.5 in ≥2 models.',fontsize=6.7)
    fig.text(.13,.114,'462 candidates are shown. Four additional candidates lack data in the other system:',fontsize=6.7)
    fig.text(.13,.085,'3 cell-line candidates + 1 organoid candidate = 466 candidates in the full union.',fontsize=6.7)
    fig.text(.13,.056,'Distinct criteria and gene-specific denominators; no joint probability or independent validation.',fontsize=6.7)
    save(fig,out/'cell_line_organoid_integration','COREAD cell-line and organoid integration')
    print('Two minimal quantitative figures generated; original figure versions remain archived.')
if __name__=='__main__':main()
