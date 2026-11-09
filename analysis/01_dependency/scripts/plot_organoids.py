"""Recreate the two quantitative SVG figures from the accompanying source CSVs.

Requires Python 3.9+, NumPy and Matplotlib. No values are jittered or normalised.
The source CSVs were checked against the final supplementary workbook.
"""
from pathlib import Path
import csv
import json
import math
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

import os
REPO = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("COREAD_FIGURE_DIR", str(REPO / "build" / "figures")))
ROOT.mkdir(parents=True, exist_ok=True)
SOURCE = REPO / "figures" / "source_data"
mpl.rcParams['svg.fonttype'] = 'none'  # Keep text editable; no font files embedded.
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.linewidth'] = 0.65
mpl.rcParams['svg.hashsalt'] = 'COREAD_1720_final_verified'


def load_rows(name):
    with (SOURCE / name).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def is_true(value):
    return value == 'True'


def save_svg(fig, name, description, dpi=160):
    svg = ROOT / (name + '.svg')
    fig.savefig(svg, format='svg', metadata={
        'Title': name.replace('_', ' '),
        'Description': description,
        'Creator': 'Matplotlib; quantitative plots of the verified COREAD supplementary table',
    })
    fig.savefig(ROOT / (name + '.png'), dpi=dpi)
    # Validate that even heatmap cells and the colour key remain vector objects.
    tree = ET.parse(svg)
    ns = {'s': 'http://www.w3.org/2000/svg'}
    assert not tree.findall('.//s:image', ns), 'Unexpected raster content in SVG'
    assert tree.findall('.//s:text', ns), 'Editable text was not preserved'
    return svg


def make_scatter():
    data = load_rows('scatter_source_data.csv')
    genes = np.array([r['Gene_symbol'] for r in data])
    x = np.array([float(r['Frequency_percent']) for r in data])
    y = np.array([float(r['Chronos_mean_3DO']) for r in data])
    candidate = np.array([is_true(r['Candidate']) for r in data])
    concordant = np.array([is_true(r['Concordant_GE_EB']) for r in data])
    assert len(data) == 1636 and candidate.sum() == 324 and concordant.sum() == 188

    fig = plt.figure(figsize=(183 / 25.4, 148 / 25.4))
    ax = fig.add_axes([0.135, 0.195, 0.715, 0.685])
    # Defaults determine all plotting colours; shape also encodes the categories.
    groups = [
        (~candidate, 'o', 0.18, 'Not selected (n = 1,312)'),
        (candidate & ~concordant, 'o', 0.82,
         'Chronos candidate; no added EB concordance (n = 136)'),
        (concordant, 'D', 0.90,
         'Chronos + re-estimated EB concordance (n = 188)'),
    ]
    for idx, (mask, marker, alpha, label) in enumerate(groups):
        ax.scatter(x[mask], y[mask], s=10, marker=marker, alpha=alpha,
                   linewidths=0, label=label, zorder=3 + idx, rasterized=False)
    ax.axvline(12.5, linewidth=0.7, linestyle=(0, (3, 3)), alpha=0.8, zorder=1)
    ax.annotate('Selection: ≥2/16', xy=(12.5, 1),
                xycoords=('data', 'axes fraction'), xytext=(3, 5),
                textcoords='offset points', ha='left', va='bottom', fontsize=6.3)
    ax.set_xlim(-3.5, 103.5)
    ax.set_ylim(-4.65, 0.65)
    ax.set_xticks([0, 12.5, 25, 50, 75, 100])
    ax.set_xticklabels(['0', '12.5', '25', '50', '75', '100'])
    ax.set_xticks(np.arange(0, 100.001, 6.25), minor=True)
    ax.set_yticks([-4, -3, -2, -1, 0])
    ax.tick_params(which='major', width=0.6, length=3, labelsize=7)
    ax.tick_params(axis='x', which='minor', width=0.4, length=1.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('COREAD 3DO models with Chronos < −0.5 (%)', fontsize=8, labelpad=8)
    ax.set_ylabel('Mean Chronos gene effect across 16 models', fontsize=8, labelpad=8)
    ax.grid(axis='y', alpha=0.16, linewidth=0.45)
    ax.set_axisbelow(True)

    positions = {g: i for i, g in enumerate(genes)}
    # Illustration labels, not an additional priority score. Points are not moved.
    right_labels = {
        'RAN': -4.36, 'SRSF3': -2.61, 'HSPD1': -2.15,
        'BCL2L1': -1.95, 'MYC': -1.74, 'UBE2I': -1.52,
        'PKM': -1.26, 'PTK2': -0.94,
    }
    for gene, y_label in right_labels.items():
        i = positions[gene]
        ax.annotate(gene, (x[i], y[i]), xytext=(106, y_label),
                    textcoords='data', ha='left', va='center', fontsize=6.6,
                    fontstyle='italic', annotation_clip=False,
                    arrowprops={'arrowstyle': '-', 'lw': 0.45, 'shrinkA': 2, 'shrinkB': 2})
    for gene, pos in {'APC': (23, -0.91), 'EGFR': (58, -1.38),
                      'HSPA5': (72, -1.06), 'KRAS': (86, -1.65)}.items():
        i = positions[gene]
        ax.annotate(gene, (x[i], y[i]), xytext=pos, textcoords='data',
                    fontsize=6.6, fontstyle='italic', ha='center', va='center',
                    arrowprops={'arrowstyle': '-', 'lw': 0.45, 'shrinkA': 2, 'shrinkB': 2})
    ax.legend(loc='lower left', bbox_to_anchor=(0.205, 0.085), frameon=False,
              fontsize=6.3, labelspacing=1.0, handletextpad=0.55,
              borderpad=0.2, scatterpoints=1)
    fig.text(0.135, 0.969, 'Recurrence and magnitude of fitness effects',
             fontsize=10, fontweight='bold', va='top')
    fig.text(0.135, 0.928, 'COREAD organoids | 16 3DO models | 1,636 evaluable genes',
             fontsize=7, va='top')
    fig.text(0.135, 0.088,
             '324 candidates: Chronos < −0.5 in ≥2 models. The gene-level mean is not a selection cutoff.',
             fontsize=6)
    fig.text(0.135, 0.062,
             'EB probabilities are re-estimated from the same effects; concordance is not independent validation.',
             fontsize=6)
    fig.text(0.135, 0.036,
             '84 input symbols without evaluable data are not plotted. Coordinates are not jittered.',
             fontsize=6)
    desc = ('1,636 evaluable genes from 1,720 input symbols. Each point has its original recurrence '
            'and mean Chronos effect in 16 COREAD 3DO models. 324 candidates satisfy Chronos < -0.5 '
            'in at least two models; 188 have additional same-model re-estimated EB probability > 0.5 '
            'in at least two models. These are exploratory criteria, not FDR-controlled discoveries. '
            'No differential COREAD selectivity is implied.')
    save_svg(fig, 'Figure_1_COREAD_recurrence_effect', desc, dpi=220)
    return fig


def make_heatmap():
    data = load_rows('heatmap_source_data.csv')
    fields = list(data[0])
    screens = [c for c in fields if c.startswith('SC-')]
    matrix = np.array([[float(row[s]) for s in screens] for row in data])
    genes = [row['Gene_symbol'] for row in data]
    counts = np.array([int(row['Positive_3DO_models']) for row in data])
    freqs = counts / 16
    concordant = np.array([is_true(row['Concordant_GE_EB']) for row in data])
    assert matrix.shape == (324, 16) and np.isfinite(matrix).all()
    assert np.array_equal((matrix < -.5).sum(axis=1), counts)
    assert concordant.sum() == 188
    n = len(genes)
    # The full-height vector canvas gives every one of the 324 genes a legible label.
    # All cells use a common, linear, original-unit scale, not row-wise Z scores.
    row_mm = 2.10
    width_mm, top_mm, bottom_mm = 210, 65, 18
    height_mm = n * row_mm + top_mm + bottom_mm
    left_mm, heatmap_width_mm = 35, 105
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    ax = fig.add_axes([left_mm / width_mm, bottom_mm / height_mm,
                       heatmap_width_mm / width_mm, n * row_mm / height_mm])
    mesh = ax.pcolormesh(np.arange(17), np.arange(n + 1), matrix,
                         shading='flat', vmin=math.floor(matrix.min()),
                         vmax=math.ceil(matrix.max()), linewidth=0,
                         edgecolors='face', antialiased=False, rasterized=False)
    ax.set_xlim(0, 16)
    ax.set_ylim(n, 0)
    ax.set_xticks(np.arange(16) + 0.5)
    ax.set_xticklabels(screens, fontsize=5.5, rotation=90, ha='left', va='center',
                        rotation_mode='anchor')
    ax.xaxis.tick_top()
    ax.tick_params(axis='x', length=0, pad=4)
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_yticklabels(genes, fontsize=5.35, fontstyle='italic')
    ax.tick_params(axis='y', length=0, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Aligned annotation columns are drawn in the same chart coordinate system.
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    ax.barh(np.arange(n) + 0.5, freqs * 0.24, left=1.055, height=0.7,
             transform=trans, clip_on=False, linewidth=0)
    for val in [0, 0.5, 1]:
        ax.text(1.055 + val * 0.24, -1.4, str(int(val * 100)),
                transform=trans, ha='center', va='bottom', fontsize=5.6, clip_on=False)
    ax.text(1.175, -5.5, 'Models (%)', transform=trans,
             ha='center', va='bottom', fontsize=6.3, fontweight='bold', clip_on=False)
    ax.text(1.47, -5.5, 'GE + EB', transform=trans,
             ha='center', va='bottom', fontsize=6.3, fontweight='bold', clip_on=False)
    ax.plot(np.full(concordant.sum(), 1.47), np.flatnonzero(concordant) + 0.5,
             linestyle='none', marker='D', markersize=2.7,
             transform=trans, clip_on=False, label='Same-model concordance')
    for boundary in [60, 106, 166]:
        ax.plot([0, 1.31], [boundary, boundary], transform=trans,
                 linestyle=(0, (2, 2)), linewidth=0.35, alpha=0.6, clip_on=False)
    ax.set_xlim(0, 16)  # Annotation columns never change the heatmap extent.
    ax.set_ylim(n, 0)

    cax = fig.add_axes([left_mm / width_mm, (height_mm - 31) / height_mm,
                        heatmap_width_mm / width_mm, 3.0 / height_mm])
    cb = fig.colorbar(mesh, cax=cax, orientation='horizontal')
    cb.set_ticks([-5, -4, -3, -2, -1, 0, 1, 2])
    cb.ax.tick_params(labelsize=6, width=0.45, length=2, pad=1.5)
    cb.outline.set_linewidth(0.5)
    if cb.solids is not None:
        cb.solids.set_rasterized(False)
    cax.set_title('Chronos gene effect | original scale; no row normalisation',
                   fontsize=6.8, loc='left', pad=4)
    fig.text(left_mm / width_mm, 1 - 5 / height_mm,
              'Fitness-effect profiles of 324 COREAD candidates',
              fontsize=10, fontweight='bold', va='top')
    fig.text(left_mm / width_mm, 1 - 11.5 / height_mm,
              '324 genes × 16 3DO models | ordered by recurrence, then mean Chronos effect',
              fontsize=7, va='top')
    fig.text(151 / width_mm, 1 - 23 / height_mm,
              'Right-hand diamond:\n≥2 models with both\nChronos < −0.5 and\nre-estimated EB P > 0.5\n(188 genes).',
              fontsize=6.1, va='top', linespacing=1.5)
    fig.text(left_mm / width_mm, 10 / height_mm,
              'All 324 gene labels are shown. Screen order follows the supplementary table; no clustering is imposed.',
              fontsize=6)
    fig.text(left_mm / width_mm, 5 / height_mm,
              'EB concordance uses the same effects and is not independent validation. Selection is exploratory, not FDR-controlled.',
              fontsize=6)
    desc = ('Heatmap of original Chronos gene-effect estimates for all 324 selected candidate genes '
            'and 16 COREAD organoid screens in 3DO format. Rows sorted by number of scores < -0.5 '
            'descending, then mean score ascending and gene symbol. Original screen order. '
            'Common linear scale [-5, 2], encompassing every value. No row normalisation, '
            'Z scoring, imputation, clipping or clustering. Right bars indicate score-call frequency; '
            'diamonds indicate at least two same-model Chronos/EB concordant calls. '
            'All 324 row labels remain editable vector text.')
    save_svg(fig, 'Figure_2_COREAD_heatmap_324x16', desc, dpi=100)
    return fig


if __name__ == '__main__':
    scatter = make_scatter()
    heatmap = make_heatmap()
    plt.close('all')
    print('Saved two SVG figures with editable text and no raster layers.')
