#!/usr/bin/env python3
"""Recompute keyed discovery statistics, frozen sets and Figure 4 source data.

Processed mode starts from patient/section measurements and upstream evidence.
It does not claim to refit the original spatial composition model. The original
matched-null vectors are retained explicitly; their undocumented Python set
iteration order prevents an exact fresh replay of historical sampling.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent


def csv_write(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('w', encoding='utf-8', newline='') as f:
        frame.to_csv(f, index=False)
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)


def json_write(obj, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp')
    with temp.open('w', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ix = np.flatnonzero(np.isfinite(p))
    order = ix[np.argsort(p[ix], kind='stable')]
    if len(order):
        adj = p[order] * len(order) / np.arange(1, len(order) + 1)
        q[order] = np.minimum(1, np.minimum.accumulate(adj[::-1])[::-1])
    return q


def flag(values):
    return values.astype(str).str.lower().isin(['true', 'yes', '1'])


def compare_frames(actual, expected, keys=('Gene',), rtol=1e-9, atol=1e-11):
    a = actual.set_index(list(keys)).sort_index()
    e = expected.set_index(list(keys)).sort_index()
    failures = []
    if not a.index.equals(e.index):
        return {'passed': False, 'failures': ['key mismatch'], 'checked_columns': 0}
    for col in e:
        if col not in a:
            failures.append({'column': col, 'reason': 'missing'})
            continue
        ev, av = e[col], a[col]
        if pd.api.types.is_bool_dtype(ev):
            equal = flag(av).to_numpy() == flag(ev).to_numpy()
        elif pd.api.types.is_numeric_dtype(ev):
            x = pd.to_numeric(av, errors='coerce').to_numpy(float)
            y = ev.to_numpy(float)
            equal = np.isclose(x, y, rtol=rtol, atol=atol, equal_nan=True)
        else:
            equal = av.fillna('').astype(str).to_numpy() == ev.fillna('').astype(str).to_numpy()
        if not np.all(equal):
            bad = np.flatnonzero(~equal)
            failures.append({'column': col, 'n_mismatch': len(bad),
                             'example_keys': [str(e.index[i]) for i in bad[:5]]})
    return {'passed': not failures, 'checked_columns': len(e.columns), 'rows': len(e), 'failures': failures}


def patient_effect_summaries(pat):
    rows = []
    specs = [('interface_effect', 'Interface', 'effect'),
             ('tumor_vs_stroma_effect', 'TumorStroma', 'effect'),
             ('distance_rho', 'Distance', 'rho')]
    for gene, d in pat.groupby('Gene', sort=False):
        row = {'Gene': gene}
        for source, prefix, suffix in specs:
            x = d[source].dropna().to_numpy(float)
            n = len(x)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                p = stats.wilcoxon(x, zero_method='wilcox', alternative='two-sided').pvalue if n >= 3 else np.nan
            row.update({f'{prefix}_N_patients': n,
                        f'{prefix}_median_{suffix}': np.median(x) if n >= 3 else np.nan,
                        f'{prefix}_mean_{suffix}': np.mean(x) if n >= 3 else np.nan,
                        f'{prefix}_positive_patient_fraction': np.mean(x > 0) if n >= 3 else np.nan,
                        f'{prefix}_Wilcoxon_P': p})
        rows.append(row)
    result = pd.DataFrame(rows)
    for _, prefix, _ in specs:
        result[f'{prefix}_BH_FDR'] = bh(result[f'{prefix}_Wilcoxon_P'])
    return result


def keyed_moran(section):
    pat = section.groupby(['patient', 'Gene'], as_index=False).Moran_I.mean()
    rows = []
    for gene, d in pat.groupby('Gene', sort=False):
        x = d.Moran_I.to_numpy(float)
        finite = x[np.isfinite(x)]
        # Preserve the archived test's NaN propagation. Omitting missing values
        # here would change the eligible hypothesis family (441 finite tests).
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p = stats.wilcoxon(x).pvalue if len(x) >= 3 else np.nan
        rows.append({'Gene': gene, 'Moran_N_patients': len(finite),
                     'Moran_median_I': np.median(finite) if len(finite) else np.nan,
                     'Moran_mean_I': np.mean(finite) if len(finite) else np.nan,
                     'Moran_positive_patient_fraction': np.mean(x > 0),
                     'Moran_Wilcoxon_P': p})
    result = pd.DataFrame(rows)
    result['Moran_BH_FDR'] = bh(result.Moran_Wilcoxon_P)
    return pat, result


def derive_hierarchy(d):
    d = d.copy()
    d['Core47'] = d.SingleCell_priority_tier.eq('A') & d.Composition_balanced_median_delta_log2CPM.gt(0) & d.Composition_balanced_BH_FDR.lt(.05)
    d['Spatial_TierA65'] = d.SingleCell_priority_tier.eq('A')
    d['Spatially_substantiated_Core47'] = d.Core47 & d.Moran_median_I.gt(0) & d.Moran_BH_FDR.lt(.05) & d.TumorStroma_positive_patient_fraction.ge(5/6)
    d['COREAD_preferential_dependency'] = d.CL_preference.eq('HIGHER_IN_COREAD')
    d['Stable_DNA_top30_80pct'] = d.DNA_top30_boot_frequency.ge(.8)
    anchor = d.COREAD_preferential_dependency | d.Stable_DNA_top30_80pct
    d['Primary_anchor_type'] = np.select(
        [d.COREAD_preferential_dependency & d.Stable_DNA_top30_80pct,
         d.COREAD_preferential_dependency, d.Stable_DNA_top30_80pct],
        ['COREAD-preferential dependency + stable top-30 patient alteration',
         'COREAD-preferential dependency', 'stable top-30 patient alteration'], default='none')
    d['Extended_anchor_5pct'] = d.COREAD_preferential_dependency | d.DNA_driver_CNA_frequency.ge(.05)
    d['Convergent9'] = d.Spatially_substantiated_Core47 & anchor
    d['Extended13_5pct'] = d.Spatially_substantiated_Core47 & d.Extended_anchor_5pct
    d['Integrated_class'] = np.select(
        [d.Convergent9, d.Extended13_5pct & ~d.Convergent9, d.Spatially_substantiated_Core47,
         d.Core47 & anchor, d.Core47, d.SingleCell_priority_tier.eq('A'), d.SingleCell_priority_tier.eq('B')],
        ['Convergent-9', 'Extended-only sensitivity', 'Spatially substantiated Core47',
         'Orthogonally anchored Core47 without spatial support', 'Other Core47',
         'Tier A composition-sensitive', 'Tier B'], default='Tier C')
    d['Spatial_category'] = np.select(
        [d.Interface_median_effect.gt(0) & d.Adjusted_interface_median_effect.gt(0),
         d.Interface_median_effect.gt(0), d.Interface_median_effect.lt(0) & d.TumorStroma_median_effect.gt(0),
         d.TumorStroma_median_effect.lt(0)],
        ['Interface-persistent', 'Interface-composition-associated', 'Tumor-interior enriched', 'Stroma-associated'],
        default='No dominant spatial category')
    d['CrossSystem_dependency_geomean'] = np.sqrt(d.CL_frequency * d.Org_frequency_3DO)
    components = {
        'CrossSystem_dependency_geomean': 'CrossSystem_dependency_geomean',
        'Tumor_epi_vs_nonEpi_median_delta_log2CPM': 'Tumor_epi_vs_nonEpi_effect',
        'Composition_balanced_median_delta_log2CPM': 'Composition_balanced_effect',
        'TumorStroma_median_effect': 'TumorStroma_median_effect',
        'Moran_median_I': 'Moran_median_I',
        'DNA_driver_CNA_frequency': 'DNA_driver_CNA_frequency'}
    for source, name in components.items():
        x = pd.to_numeric(d[source], errors='coerce')
        d['z_' + name] = (x - x.mean()) / x.std(ddof=1)
    d['Integrated_visual_score'] = d[['z_' + c for c in components.values()]].mean(axis=1)
    # Stable ordering is documented; no rank feeds back into candidate selection.
    d = d.sort_values('Integrated_visual_score', ascending=False, kind='stable').reset_index(drop=True)
    d['Figure_visual_rank'] = np.arange(1, len(d)+1)
    d['Figure_A_class'] = np.select([d.Convergent9, d.Core47], ['Convergent9', 'Core47'], default='Other candidate')
    d['Figure_C_Top24'] = d.Figure_visual_rank.le(24)
    d['Figure_C_row_order'] = np.where(d.Figure_C_Top24, d.Figure_visual_rank, np.nan)
    d['Figure_D_Top25'] = d.Figure_visual_rank.le(25)
    d['Figure_D_visual_rank'] = np.where(d.Figure_D_Top25, d.Figure_visual_rank, np.nan)
    heat = {
        'dependency': 'CrossSystem_dependency_geomean',
        'tumor_vs_normal': 'Tumor_vs_normal_epi_median_delta_log2CPM',
        'epi_vs_non_epi': 'Tumor_epi_vs_nonEpi_median_delta_log2CPM',
        'composition_balanced': 'Composition_balanced_median_delta_log2CPM',
        'tumor_stroma': 'TumorStroma_median_effect', 'Moran': 'Moran_median_I', 'DNA': 'DNA_driver_CNA_frequency'}
    for name, source in heat.items():
        x = d.loc[d.Figure_C_Top24, source]
        d[f'Figure_C_scaled_{name}_top24'] = np.nan
        d.loc[d.Figure_C_Top24, f'Figure_C_scaled_{name}_top24'] = (x-x.mean())/x.std(ddof=1)
    return d


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['processed','counts','raw'], default='processed')
    parser.add_argument('--raw-dir', type=Path, help='Directory containing sp.h5ad, W_cell_density_q05.csv and pathology/')
    parser.add_argument('--output-dir', type=Path, default=ROOT/'results')
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    upstream = pd.read_csv(ROOT/'data/upstream_gene_evidence.csv')
    pat = pd.read_csv(ROOT/'data/patient_gene_effects.csv')
    section = pd.read_csv(ROOT/'data/moran_section_466.csv')
    source_status = None
    if args.mode in ['counts','raw']:
        source_out = out/'source_reconstruction'
        command = [sys.executable, str(ROOT/'scripts/recompute_from_raw.py'),
                   '--upstream-evidence', str(ROOT/'data/upstream_gene_evidence.csv'),
                   '--output-dir', str(source_out)]
        if args.mode == 'counts':
            command += ['--checkpoint-dir', str(ROOT/'data/counts')]
        else:
            if args.raw_dir is None: parser.error('--raw-dir is required in raw mode')
            command += ['--h5ad', str(args.raw_dir/'sp.h5ad'),
                        '--annotations-dir', str(args.raw_dir/'pathology'),
                        '--cell-density', str(args.raw_dir/'W_cell_density_q05.csv')]
        subprocess.run(command, check=True)
        source_status = json.loads((source_out/'status.json').read_text())
        source_table = pd.read_csv(source_out/'ST11_from_source_counts.csv').set_index('Gene')
        refreshed = [c for c in upstream if c.startswith(('Adjusted_interface_','Spatial_mean_','Spatial_detection_')) or c=='Spatial_module']
        upstream = upstream.set_index('Gene')
        upstream[refreshed] = source_table.reindex(upstream.index)[refreshed]
        upstream = upstream.reset_index()
        pat = pd.read_csv(source_out/'patient_gene_effects.csv')
        section = pd.read_csv(source_out/'moran_section_466.csv')
    expected = pd.read_csv(ROOT/'reference/discovery_spatial_expected.csv')
    assert upstream.Gene.nunique() == 466 and len(upstream) == 466
    effects = patient_effect_summaries(pat)
    moran_pat, moran = keyed_moran(section)
    combined = upstream.merge(effects, on='Gene', validate='one_to_one').merge(moran, on='Gene', validate='one_to_one')
    # Matching bins preserve original catalogue order before visual sorting.
    expr = pd.qcut(combined.Spatial_mean_logCPM.rank(method='first'), 5, labels=False)
    det = pd.qcut(combined.Spatial_detection_fraction.rank(method='first'), 5, labels=False)
    combined['match_bin'] = expr.astype(str) + '_' + det.astype(str)
    result = derive_hierarchy(combined)
    result = result.merge(pd.read_csv(ROOT/'config/figure_label_selection.csv'), on='Gene', validate='one_to_one')
    result = result[expected.columns]
    verification = compare_frames(result, expected)
    csv_write(result, out/'discovery_spatial_recomputed.csv')
    csv_write(moran_pat, out/'patient_Moran_keyed.csv')
    csv_write(moran, out/'Moran_gene_tests_keyed.csv')
    csv_write(effects, out/'spatial_patient_effect_summaries.csv')
    sets = {'Core47': result.loc[result.Core47, 'Gene'].tolist(),
            'Spatial43': result.loc[result.Spatially_substantiated_Core47, 'Gene'].tolist(),
            'Convergent9': result.loc[result.Convergent9, 'Gene'].tolist(),
            'Extended13': result.loc[result.Extended13_5pct, 'Gene'].tolist()}
    json_write(sets, out/'frozen_gene_sets.json')
    null = np.load(ROOT/'data/core47_matched_permutations_nulls.npz')
    rows = []
    for metric in null.files:
        arr = null[metric]
        observed = float(result.loc[result.Core47, metric].mean())
        exceed = int(np.sum(arr >= observed))
        rows.append({'metric': metric, 'observed_core47_mean': observed, 'null_mean': float(arr.mean()),
                     'null_sd': float(arr.std(ddof=1)), 'null_p2.5': float(np.quantile(arr,.025)),
                     'null_p97.5': float(np.quantile(arr,.975)), 'empirical_P_greater': (exceed+1)/(len(arr)+1),
                     'exceedances': exceed})
    null_summary = pd.DataFrame(rows)
    csv_write(null_summary, out/'core47_matched_permutations_replayed.csv')
    null_check = compare_frames(null_summary, pd.read_csv(ROOT/'reference/core47_matched_permutations_summary.csv'), keys=('metric',))
    sources = out/'source_data'
    csv_write(result[['Gene', 'Figure_A_class','Figure_A_label','Extended13_5pct','CrossSystem_dependency_geomean',
                      'Composition_balanced_median_delta_log2CPM','DNA_driver_CNA_frequency']], sources/'Figure_4A.csv')
    csv_write(result.loc[result.Core47, ['Gene','Spatial_category','Convergent9','Figure_B_label',
                      'TumorStroma_median_effect','Moran_median_I','CrossSystem_dependency_geomean']], sources/'Figure_4B.csv')
    csv_write(result.loc[result.Figure_C_Top24, ['Gene','Figure_C_row_order'] + [c for c in result if c.startswith('Figure_C_scaled_')] +
                         ['Immune_axis','Metabolic_axis','Endocrine_axis','Neural_axis','Stromal_matrix_axis']], sources/'Figure_4C.csv')
    csv_write(result.loc[result.Figure_D_Top25, ['Gene','Integrated_visual_score','Integrated_class','Figure_D_visual_rank']], sources/'Figure_4D.csv')
    counts = {name: len(genes) for name, genes in sets.items()}
    required = {'Core47': 47, 'Spatial43': 43, 'Convergent9': 9, 'Extended13': 13}
    numerical_checks = {
        'hierarchy_counts': counts == required,
        'Moran_finite_hypotheses_441': int(result.Moran_Wilcoxon_P.notna().sum()) == 441,
        'unique_genes_466': result.Gene.nunique() == 466,
        'C9_members': set(sets['Convergent9']) == {'BCL2L1','CTNNB1','FASN','HMGB1','IRS2','MYC','RNF6','SOX9','YAP1'},
        'Core47_Moran_positive_significant_44': int((result.Core47 & result.Moran_median_I.gt(0) & result.Moran_BH_FDR.lt(.05)).sum()) == 44,
        'sets_nested': set(sets['Convergent9']) <= set(sets['Extended13']) <= set(sets['Spatial43']) <= set(sets['Core47'])}
    if not args.no_plots:
        from plot_figure4 import plot_figure4
        plot_figure4(result, out/'figures')
    status = {'module': '04_discovery_spatial', 'mode': args.mode,
              'passed': bool(verification['passed'] and null_check['passed'] and all(numerical_checks.values())),
              'counts': counts, 'reference_comparison': verification,
              'matched_null_replay_comparison': null_check, 'numerical_checks': numerical_checks,
              'execution_scope': 'Patient/section checkpoints to keyed inference, hierarchy and plots' if args.mode=='processed' else 'Source counts through regional measurements, composition refit, Moran, hierarchy and plots',
              'source_reconstruction': source_status,
              'limits': (['Composition-adjusted interface summaries and spatial coexpression modules are retained source checkpoints; use counts/raw mode to refit them.'] if args.mode=='processed' else ['Source pathology and cell2location factors are reused; source deconvolution and original sequencing alignment are not refitted.']) +
                        ['Reported discovery matched-null statistics replay archived vectors; exact historical new-draw sampling cannot be reconstructed because Python set iteration order was not saved.']}
    json_write(status, out/'status.json')
    print(json.dumps(status, indent=2))
    if not status['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
