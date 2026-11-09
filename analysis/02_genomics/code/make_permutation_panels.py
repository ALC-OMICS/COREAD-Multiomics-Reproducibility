"""Compact, source-exact summary of the 15 previously displayed concordance tests.
The primary driver/CNA-union tests and arm sensitivities remain in Table 9 and
original Figures 5/6. No FDR family is recalculated on this visual subset.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures/permutation_summary';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'svg.fonttype':'none','font.family':'DejaVu Sans'})
with (ROOT/'results/permutation_statistics.csv').open(newline='') as f:stats=list(csv.DictReader(f))
cohorts=['AC_ICAM','CPTAC2','PANCAN']
labels={'Any_exported_mutation':'Coding/splice mutations','Driver_mutation':'Putative-driver mutations',
        'Amplification':'Amplification','Deep_deletion':'Deep deletion','CNA_any':'Amplification or deep deletion',
        'DNA_broad':'Any mutation or high-level CNA'}
rows=[]
for co in cohorts:
    eps=['Any_exported_mutation' if co=='AC_ICAM' else 'Driver_mutation',
         'Amplification','Deep_deletion','CNA_any','DNA_broad']
    for ep in eps:
        found=[r for r in stats if r['Reference']=='COREAD_GDC' and r['Comparison']==co and
               r['Endpoint']==ep and r['Null_scheme']=='GENE_LABEL_PERMUTATION']
        if len(found)!=1 or not found[0]['Empirical_P_greater']:raise ValueError((co,ep))
        rows.append(found[0])
fig=plt.figure(figsize=(8.7,7.5))
ax=fig.add_axes([.33,.18,.425,.68])
y=np.array([j+(j//5)*1.25 for j in range(15)])
for pos,r in zip(y,rows):
    ax.hlines(pos,float(r['Null_percentile_2p5']),float(r['Null_percentile_97p5']),
              linewidth=3,alpha=.22)
ax.scatter([float(r['Spearman_rho']) for r in rows],y,s=25,label='Observed Spearman rho',zorder=3)
ax.axvline(0,linestyle=':',linewidth=.8,alpha=.6)
ax.set_yticks(y);ax.set_yticklabels([labels[r['Endpoint']] for r in rows],fontsize=8)
ax.set(xlim=(-.19,1.04),ylim=(max(y)+.65,-1.1),xticks=[0,.25,.5,.75,1])
ax.tick_params(axis='y',length=0,pad=7);ax.tick_params(axis='x',labelsize=8.5)
ax.spines[['top','right','left']].set_visible(False)
ax.set_xlabel('Spearman rank concordance',fontsize=10,labelpad=8)
for idx,co in enumerate(cohorts):
    yy=y[idx*5]-.8
    h={'AC_ICAM':'AC-ICAM | external colon cohort','CPTAC2':'CPTAC-2 | external colon cohort',
       'PANCAN':'PanCancer | shared-TCGA robustness'}[co]
    ax.text(-.68,yy,h,transform=ax.get_yaxis_transform(),fontweight='bold',fontsize=8.5,
            va='center',clip_on=False)
    if idx:ax.axhline(yy-.5,linewidth=.6,alpha=.18)
for pos,r in zip(y,rows):
    ax.text(1.115,pos,f"{float(r['Empirical_P_greater']):.1e}",transform=ax.get_yaxis_transform(),
            ha='left',va='center',fontsize=8)
    ax.text(1.405,pos,f"{float(r['BH_FDR']):.1e}",transform=ax.get_yaxis_transform(),
            ha='left',va='center',fontsize=8)
ax.text(1.115,-1.62,'Empirical P',transform=ax.get_yaxis_transform(),fontsize=8.3,fontweight='bold')
ax.text(1.405,-1.62,'BH FDR',transform=ax.get_yaxis_transform(),fontsize=8.3,fontweight='bold')
fig.text(.04,.96,'Permutation concordance across cohorts',fontsize=13,fontweight='bold')
fig.text(.04,.927,'10,000 permutations per test; all matched genes, not only the displayed top 30.',fontsize=9)
fig.legend(handles=[Line2D([],[],marker='o',linestyle='none',markersize=5,label='Observed concordance'),
                    Line2D([],[],linewidth=3,alpha=.22,label='Central 95% permutation null range')],
           loc='lower left',bbox_to_anchor=(.035,.086),frameon=False,ncol=2,fontsize=8.3)
fig.text(.04,.069,'Null ranges are not confidence intervals for observed correlations. Original BH families are retained.',fontsize=7.8)
fig.text(.04,.046,'AC-ICAM has no driver annotations. Broad DNA composition is distinct from the primary driver/CNA ranking.',fontsize=7.6)
fig.text(.04,.023,'This compact selection omits the driver/CNA-union and arm-sensitivity tests; complete results: Table 9.',fontsize=7.6)
stem='Figure_8_Permutation_summary_compact'
fig.savefig(OUT/(stem+'.svg'));fig.savefig(OUT/(stem+'.png'),dpi=180);plt.close(fig)
with (OUT/(stem+'_source.csv')).open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('SAVED',stem,'(15 exact source contrasts; all empirical P and BH values shown)')
