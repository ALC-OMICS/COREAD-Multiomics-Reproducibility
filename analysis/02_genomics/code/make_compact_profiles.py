"""Rebuild compact composition plots from the exact CSV counterparts of Tables 7/8.

This is a reconstructed plotting script, not a new statistical analysis. Bar widths
are unit-sum normalized marginal frequencies, NOT a patient-level alteration union.
The right-hand DNA percentages are separate, measured composite frequencies.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'figures/compact_profiles'
OUT.mkdir(parents=True, exist_ok=True)
# Colours requested through the user's reference graphic.
COLORS = {'mutation':'#c9c86b', 'amp':'#d3a33a', 'del':'#f28e2b',
          'rna_high':'#d26b6b', 'rna_low':'#5a8fde'}
CATS = list(COLORS)
plt.rcParams.update({'svg.fonttype':'none', 'font.family':'DejaVu Sans'})


def read(name):
    with (ROOT/'results'/name).open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))


def records(cohort):
    primary=read('primary_466_results.csv')
    top=[r['Gene'] for r in sorted(primary,key=lambda r:int(r['DNA_display_rank']))[:30]]
    if cohort=='TCGA':
        source={r['Gene']:r for r in primary}
        cols=['Driver_mutation_frequency','Amplification_frequency','Deep_deletion_frequency',
              'RNA_high_frequency','RNA_low_frequency']
        nn=['Driver_mutation_n_profiled','Amplification_n_profiled','Deep_deletion_n_profiled',
            'RNA_high_n_profiled','RNA_low_n_profiled']
        f,k,n='DNA_driver_CNA_frequency','DNA_driver_CNA_altered','DNA_driver_CNA_n_profiled'
    else:
        source={r['Gene']:r for r in read('external_466_results_wide.csv') if r['Cohort']==cohort}
        mutation='Any_exported_mutation' if cohort=='AC_ICAM' else 'Driver_mutation'
        cols=[mutation+'_Frequency','Amplification_Frequency','Deep_deletion_Frequency',
              'RNA_high_Frequency','RNA_low_Frequency']
        nn=[mutation+'_N','Amplification_N','Deep_deletion_N','RNA_high_N','RNA_low_N']
        composite='DNA_broad' if cohort=='AC_ICAM' else 'DNA_driver_CNA'
        f,k,n=composite+'_Frequency',composite+'_Altered',composite+'_N'
    rows=[]
    for rank,g in enumerate(top,1):
        r=source[g]
        if any(not r[c] for c in cols+[f,k,n]):
            raise ValueError(f'{cohort}:{g}: missing data; do not silently replace with zero.')
        values=np.array([float(r[c]) for c in cols])
        if np.any((values<0)|(values>1)) or values.sum()<=0:
            raise ValueError(f'{cohort}:{g}: invalid or undefined composition')
        row={'Cohort':cohort,'TCGA_display_rank':rank,'Gene':g,'DNA_endpoint':f.rsplit('_',1)[0],
             'DNA_n':int(float(r[n])),'DNA_altered':int(float(r[k])),
             'DNA_frequency':float(r[f]),'Normalization_sum_of_marginal_frequencies':float(values.sum())}
        for j,cat in enumerate(CATS):
            row[cat+'_n']=int(float(r[nn[j]]))
            row[cat+'_frequency']=float(values[j])
            row[cat+'_composition']=float(values[j]/values.sum())
            row[cat+'_observed']=bool(values[j]>0)
        rows.append(row)
    return rows


def render(cohort):
    rows=records(cohort)
    name={'TCGA':'TCGA COREAD','CPTAC2':'CPTAC-2','AC_ICAM':'AC-ICAM'}[cohort]
    fig=plt.figure(figsize=(8.0,10.8))
    ax=fig.add_axes([.18,.155,.79,.744])
    y=np.arange(30);left=np.zeros(30);h=.72
    for cat in CATS:
        widths=np.array([r[cat+'_composition'] for r in rows])
        ax.barh(y,widths,left=left,height=h,color=COLORS[cat],edgecolor='white',linewidth=.25)
        left+=widths
    for i,r in enumerate(rows):
        for j,cat in enumerate(CATS):
            ax.add_patch(Rectangle((-.143+j*.023,i-h/2),.019,h,
                facecolor=COLORS[cat] if r[cat+'_observed'] else 'white',edgecolor='.55',linewidth=.45))
        ax.text(1.132,i,f"{100*r['DNA_frequency']:.1f}%",ha='right',va='center',fontsize=9)
    np.testing.assert_allclose(left,1,atol=1e-14)
    ax.set(xlim=(-.165,1.16),ylim=(29.65,-1.0),yticks=y,xticks=np.linspace(0,1,6))
    ax.set_yticklabels([r['Gene'] for r in rows],fontstyle='italic',fontsize=9)
    ax.tick_params(axis='y',length=0,pad=6)
    ax.tick_params(axis='x',labelsize=8.6,length=3)
    ax.spines[['top','right','left']].set_visible(False)
    ax.spines['bottom'].set_bounds(0,1)
    ax.axvline(0,linewidth=.65,color='.3')
    ax.set_xlabel('Normalized alteration-type composition',fontsize=10,labelpad=7)
    ax.text(-.087,-.83,'Types',ha='center',va='bottom',fontsize=7.5)
    ax.text(1.132,-.83,'DNA altered*',ha='right',va='bottom',fontsize=8.2)
    fig.text(.065,.965,name,fontsize=15,fontweight='bold')
    mutation='Coding/splice mutation' if cohort=='AC_ICAM' else 'Putative-driver mutation'
    fig.text(.065,.935,'Same 30 genes and order, fixed by the TCGA primary DNA ranking.',fontsize=8.5)
    if cohort=='AC_ICAM':
        subtitle='Mutation/CNA n = 281; RNA n = 348. Driver classifications unavailable.'
        footer='*Any coding/splice mutation OR amplification OR deep deletion; 281 DNA-profiled patients.'
    elif cohort=='CPTAC2':
        subtitle='Mutation n = 106; CNA n = 105; RNA n = 106; joint DNA profiles n = 105.'
        footer='*Putative-driver mutation OR amplification OR deep deletion; 105 joint-profiled patients.'
    else:
        subtitle='Mutation n = 581; CNA n = 621; RNA n = 624; joint DNA profiles n = 572.'
        footer='*Putative-driver mutation OR amplification OR deep deletion; 572 joint-profiled patients.'
    fig.text(.065,.915,subtitle,fontsize=8)
    labels=[mutation,'CNA amplification','CNA deep deletion','mRNA high','mRNA low']
    handles=[Rectangle((0,0),1,1,facecolor=COLORS[c],edgecolor='none') for c in CATS]
    fig.legend(handles,labels,loc='lower left',bbox_to_anchor=(.057,.071),ncol=3,
               frameon=False,fontsize=8.0,columnspacing=1.6,handlelength=1.1,handletextpad=.45)
    fig.text(.065,.058,footer,fontsize=7.5)
    fig.text(.065,.038,'Bars show relative composition, not patient prevalence; right labels are separate DNA frequencies.',fontsize=7.3)
    fig.text(.065,.019,'Mini-tiles: alteration-type presence, not biological axes. RNA is source-specific; protein/SV not evaluated.',fontsize=7.2)
    stem=f'Figure_{cohort}_compact_alteration_profile_top30'
    fig.savefig(OUT/(stem+'.svg'))
    fig.savefig(OUT/(stem+'.png'),dpi=180)
    plt.close(fig)
    with (OUT/(stem+'_source.csv')).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('SAVED',stem)
    return stem


def main():
    stems=[render(c) for c in ['TCGA','CPTAC2','AC_ICAM']]
    (OUT/'README.md').write_text('''# Compact alteration-type composition profiles

These numerical plots reproduce the later requested stacked-bar representation,
not the original patient-by-gene OncoPrint. Source values are the exact CSV
counterparts of Supplementary Tables 7 and 8. The same TCGA top 30 is retained.

For gene g and displayed modality m, bar width is f(g,m) / sum_m f(g,m), using
five marginal frequencies: mutation, amplification, deep deletion, RNA high and
RNA low. They use their original modality-specific denominators. All five
frequencies are finite for these 30 genes in each cohort. The sum is a descriptive
normalization; it is NOT a proportion of unique altered patients, mutation burden,
mean frequency or validation statistic. RNA thresholds are not harmonized.

Right percentages remain unnormalized unique-patient DNA-union frequencies:
TCGA/CPTAC use putative-driver OR amplification OR deep deletion; AC-ICAM uses
broad coding/splice mutation OR amplification OR deep deletion, not drivers.
Thus the right-hand endpoint is not identically defined in AC-ICAM. Mini-tiles
indicate observed modality calls and do not encode membership of the five axes.

No protein or structural-variant segments are fabricated. The layout and labels
were clarified for this repository; the delivered numeric SVGs are preserved
under archive/delivered_figure_versions/compact_profiles/. AI-generated versions
are excluded because they do not faithfully represent numerical source values.
''',encoding='utf-8')


if __name__=='__main__':main()
