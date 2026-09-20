# -*- coding: utf-8 -*-
"""
regenerate_fig4C.py
Regenerates panel C of Figure4_DEG_pathway, and the complete Figure 4, fixing the
defects raised by the Chinese team on 2026-08-12 plus two they did not report.

Defects fixed
-------------
(1) [reported]  Pathway names not fully displayed.
    -> every panel now carries the readable name AND the full MSigDB Hallmark
       identifier, on two lines, with no clipping at any output size.

(2) [reported]  In the 'Inflammatory response' panel the 0 gridline did not sit at
    the foot of the bars; the axis started below 0.
    -> cause: SCART1 has a genuinely negative group mean (-0.00096).  matplotlib
       auto-extended ylim below zero, so the 0 tick floated above the axis floor.
    -> fix: a single shared y-axis with an explicit, hand-set tick sequence and a
       solid zero baseline that the bars visibly stand on.  The negative value is
       still drawn (it is real data) but is now unmistakably read against zero.

(3) [not reported] Duplicated y-tick labels.  In the submitted panel the
    'IFN-gamma response' axis printed 0.12 twice and 'IL6/JAK/STAT3' printed 0.14
    twice, because five independent auto-scaled axes were rendered at 7 pt into a
    narrow gridspec cell.  The shared axis with explicit ticks removes this class
    of bug entirely.

(4) [not reported] Five independent y-scales made the five pathways look directly
    comparable when they were not.  A shared scale makes the comparison honest.

Two versions are produced:
    Figure4C_v1_bars        - faithful drop-in replacement (bar = mean, whisker = SEM)
    Figure4C_v2_distribution- recommended; adds the per-spot distribution behind
                              each bar, which pre-empts the reviewer objection that
                              a SEM computed over spots from n = 1 mouse per arm is
                              not an inferential error bar.

Outputs PNG (600 dpi), PDF (vector, editable text) and TIFF (600 dpi, LZW).
"""
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm

# ---------------------------------------------------------------- paths ----
DATA    = Path(os.environ.get('SCART_DATA', './data'))
RESULTS = Path(os.environ.get('RESULTS', './results'))
QC = RESULTS / 'qc'
OUT = RESULTS / 'figures'
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

GROUPS = ['control', 'SBRT', 'SCART1', 'SCART3']
GROUP_LABEL = {'control': 'Control', 'SBRT': 'SBRT',
               'SCART1': 'SCART\u00d71', 'SCART3': 'SCART\u00d73'}
GC = {'control': '#7f7f7f', 'SBRT': '#4575b4', 'SCART1': '#fdae61', 'SCART3': '#d73027'}

# MSigDB id, readable name (pre-wrapped so nothing can clip or collide), id suffix
PATHWAYS = [
    ('HALLMARK_INTERFERON_GAMMA_RESPONSE', 'Interferon-\u03b3\nresponse',       'INTERFERON_GAMMA\n_RESPONSE'),
    ('HALLMARK_INFLAMMATORY_RESPONSE',     'Inflammatory\nresponse',            'INFLAMMATORY\n_RESPONSE'),
    ('HALLMARK_TNFA_SIGNALING_VIA_NFKB',   'TNF\u03b1 signalling\nvia NF-\u03baB', 'TNFA_SIGNALING\n_VIA_NFKB'),
    ('HALLMARK_IL6_JAK_STAT3_SIGNALING',   'IL-6 / JAK / STAT3\nsignalling',    'IL6_JAK_STAT3\n_SIGNALING'),
    ('HALLMARK_COMPLEMENT',                'Complement\n',                      'COMPLEMENT\n'),
]

# ---------------------------------------------------------------- data -----
pw = pd.read_csv(QC / 'pathway_aucell_per_cell.csv', index_col=0)
present = [g for g in GROUPS if g in set(pw['phase2_group'])]
N = pw['phase2_group'].value_counts()

stat = {}
for col, _, _ in PATHWAYS:
    g = pw.groupby('phase2_group')[col]
    stat[col] = pd.DataFrame({'mean': g.mean(), 'sem': g.sem(),
                              'median': g.median(), 'n': g.size()}).reindex(present)

# One shared scale for all five panels.  Two variants: the bar version only has to
# hold mean+SEM; the distribution version has to hold the plotted spread as well,
# otherwise the violins are silently guillotined at the axis top.
_bar_top = max((stat[c]['mean'] + stat[c]['sem']).max() for c, _, _ in PATHWAYS)
PCT_LO, PCT_HI = 2.5, 97.5      # violins are trimmed to this range, and said so
_dist_top = max(np.percentile(pw[c].dropna(), PCT_HI) for c, _, _ in PATHWAYS)

YMIN = -0.02
YMAX_BAR = np.ceil(_bar_top * 20) / 20 + 0.01
YMAX_DIST = np.ceil(_dist_top * 10) / 10
print(f'  y-max: bars {YMAX_BAR:.2f}, distribution {YMAX_DIST:.2f} '
      f'({PCT_HI:.1f}th pct = {_dist_top:.3f})')


def style_axis(ax, k, ymax, step):
    """Shared, explicit, non-duplicating y-axis with a real zero baseline."""
    ticks = np.round(np.arange(0.0, ymax + 1e-9, step), 2)
    ax.set_ylim(YMIN, ymax)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{v:.2f}' for v in ticks] if k == 0 else [])
    ax.spines['bottom'].set_visible(False)          # axis floor is not zero
    ax.axhline(0, color='black', lw=0.9, zorder=4)  # zero IS the baseline
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([GROUP_LABEL[g] for g in present],
                       rotation=40, ha='right', fontsize=6.5)
    ax.tick_params(axis='x', length=0, pad=1)


def panel_title(ax, name, ident):
    """Readable name above, MSigDB suffix below, both pre-wrapped so that neither
    clips at the panel edge nor collides with the neighbouring panel."""
    ax.set_title(name, fontsize=7, pad=13, linespacing=1.35)
    ax.text(0.5, 1.015, ident, transform=ax.transAxes, ha='center', va='bottom',
            fontsize=4.9, color='#666666', linespacing=1.25)


def draw_C(axes, distribution=False):
    ymax = YMAX_DIST if distribution else YMAX_BAR
    step = 0.10 if distribution else 0.05
    for k, (col, name, ident) in enumerate(PATHWAYS):
        ax = axes[k]
        s = stat[col]
        if distribution:
            data = []
            for g in present:
                d = pw.loc[pw['phase2_group'] == g, col].dropna().values
                lo, hi = np.percentile(d, [PCT_LO, PCT_HI])
                data.append(d[(d >= lo) & (d <= hi)])
            parts = ax.violinplot(data, positions=range(len(present)), widths=0.82,
                                  showextrema=False, showmedians=False)
            for b, g in zip(parts['bodies'], present):
                b.set_facecolor(GC[g]); b.set_alpha(0.22)
                b.set_edgecolor(GC[g]); b.set_linewidth(0.4)
            ax.bar(range(len(present)), s['mean'].values, width=0.42,
                   color=[GC[g] for g in present], edgecolor='black',
                   linewidth=0.45, zorder=3)
            for i, g in enumerate(present):
                ax.plot([i - 0.21, i + 0.21], [s['median'].iloc[i]] * 2,
                        color='black', lw=0.9, zorder=5)
        else:
            ax.bar(range(len(present)), s['mean'].values,
                   yerr=s['sem'].values, width=0.68,
                   color=[GC[g] for g in present], edgecolor='black',
                   linewidth=0.5, capsize=2, zorder=3,
                   error_kw=dict(elinewidth=0.6, capthick=0.6, zorder=5))
        style_axis(ax, k, ymax, step)
        panel_title(ax, name, ident)
        if k == 0:
            ax.set_ylabel('Hallmark programme score\n(Scanpy score_genes, per spot)',
                          fontsize=6.8)


def footnote(fig, y, distribution):
    n_txt = ',  '.join(f'{GROUP_LABEL[g]} n={N[g]}' for g in present)
    if distribution:
        msg = (f'Bar = group mean;  horizontal rule = median;  shaded silhouette = per-spot '
               f'distribution trimmed to the {PCT_LO:g}–{PCT_HI:g}th percentile '
               f'(a small number of high-scoring spots therefore fall outside the plotted '
               f'range).  Rim (L-zone) PTPRC+ spots only:  ' + n_txt + '.')
    else:
        msg = ('Bar = group mean;  whisker = SEM across spots.  Rim (L-zone) PTPRC+ spots '
               'only:  ' + n_txt + '.')
    msg += ('  Grey text under each title is the MSigDB Hallmark identifier '
            '(all prefixed HALLMARK_).')
    msg += ('\nScores are Scanpy score_genes values (expression relative to a matched '
            'background gene set), so zero is meaningful and negative values are possible; '
            'the SCART\u00d71 inflammatory-response mean is \u22120.001.'
            '\nn = 1 mouse per arm \u2014 spots are not independent replicates, so these '
            'comparisons are descriptive and no inferential test is reported.')
    fig.text(0.005, y, msg, fontsize=5.6, va='top', linespacing=1.45, color='#333333')


def save(fig, stem):
    for ext, kw in [('png', dict(dpi=600)), ('pdf', {}),
                    ('tiff', dict(dpi=600, pil_kwargs={'compression': 'tiff_lzw'}))]:
        fig.savefig(OUT / f'{stem}.{ext}', facecolor='white',
                    bbox_inches='tight', pad_inches=0.04, **kw)
    plt.close(fig)
    print(f'  [saved] {stem}.{{png,pdf,tiff}}')


# =============== standalone panel C, two versions ==========================
print('[1] standalone panel C')
for stem, dist in [('Figure4C_v1_bars', False), ('Figure4C_v2_distribution', True)]:
    fig, axes = plt.subplots(1, 5, figsize=(7.2, 3.05))
    draw_C(axes, distribution=dist)
    plt.tight_layout(rect=[0, 0.20, 1, 0.965])
    footnote(fig, 0.175, dist)
    save(fig, stem)

# =============== complete Figure 4 with panel C fixed ======================
print('[2] complete Figure 4 (A and B unchanged, C rebuilt)')
deg = pd.read_csv(QC / 'deg_SCART3L_vs_SBRTL.csv')
panel_means = pd.read_csv(QC / 'phase1_panel_means_in_L_zone.csv', index_col=0)

fig = plt.figure(figsize=(7.2, 7.8))
gs = gridspec.GridSpec(3, 5, figure=fig, height_ratios=[1.15, 1.15, 1.0],
                       hspace=0.75, wspace=0.30,
                       left=0.10, right=0.955, top=0.930, bottom=0.185)

# --- 4A volcano -----------------------------------------------------------
ax = fig.add_subplot(gs[0, 0:2])
ax.text(-0.22, 1.06, 'A', transform=ax.transAxes, fontsize=10, fontweight='bold', va='top')
v = deg.copy()
v['mlp'] = -np.log10(v['pvals_adj'].clip(lower=1e-300))
v.loc[v['mlp'] > 50, 'mlp'] = 50
sig = (v['pvals_adj'] < 0.05) & (v['logfoldchanges'].abs() > 0.5)
ax.scatter(v.loc[~sig, 'logfoldchanges'], v.loc[~sig, 'mlp'], s=2.5,
           c='#cccccc', alpha=0.45, rasterized=True, linewidths=0)
ax.scatter(v.loc[sig & (v.logfoldchanges > 0), 'logfoldchanges'],
           v.loc[sig & (v.logfoldchanges > 0), 'mlp'], s=5, c='#d73027',
           alpha=0.85, label=f'\u2191 SCART\u00d73 (n={int(N["SCART3"])})',
           rasterized=True, linewidths=0)
ax.scatter(v.loc[sig & (v.logfoldchanges < 0), 'logfoldchanges'],
           v.loc[sig & (v.logfoldchanges < 0), 'mlp'], s=5, c='#4575b4',
           alpha=0.85, label=f'\u2191 SBRT (n={int(N["SBRT"])})',
           rasterized=True, linewidths=0)
lab_rows = pd.concat([v.loc[sig & (v.logfoldchanges > 0)].nlargest(8, 'logfoldchanges'),
                      v.loc[sig & (v.logfoldchanges < 0)].nsmallest(8, 'logfoldchanges')])
texts = [ax.text(r.logfoldchanges, r.mlp, r.names, fontsize=5)
         for r in lab_rows.itertuples()]
try:                                   # de-collide the gene labels (defect in v1)
    from adjustText import adjust_text
    adjust_text(texts, ax=ax, expand=(1.25, 1.5),
                arrowprops=dict(arrowstyle='-', color='#888888', lw=0.3))
except Exception as e:
    print('   [warn] adjustText unavailable, labels left as-is:', e)
ax.axhline(-np.log10(0.05), color='grey', ls='--', lw=0.4)
ax.axvline(0.5, color='grey', ls=':', lw=0.4)
ax.axvline(-0.5, color='grey', ls=':', lw=0.4)
ax.set_xlabel('log\u2082 fold-change (SCART\u00d73 / SBRT)')
ax.set_ylabel('\u2212log\u2081\u2080 (adjusted p)')
ax.set_title('DEG \u2014 SCART\u00d73 Rim vs SBRT Rim (PTPRC+ immune)', fontsize=7, pad=4)
ax.legend(fontsize=6, loc='upper left', frameon=False)

# --- 4B immune panel heatmap ---------------------------------------------
ax = fig.add_subplot(gs[0, 2:5])
ax.text(-0.12, 1.06, 'B', transform=ax.transAxes, fontsize=10, fontweight='bold', va='top')
pm = panel_means.copy()
pm.columns = [c.replace('_L', '') for c in pm.columns]
pm = pm.reindex(columns=GROUPS)
z = pm.sub(pm.mean(axis=1), axis=0).div(pm.std(axis=1).replace(0, 1), axis=0)
im = ax.imshow(z.values, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
ax.set_xticks(range(z.shape[1]))
ax.set_xticklabels([GROUP_LABEL[c] for c in z.columns], fontsize=6.5)
ax.set_yticks(range(z.shape[0]))
ax.set_yticklabels(z.index, fontsize=6)
ax.tick_params(axis='y', pad=1, length=2)
ax.set_title('Phase-1 immune gene panel \u2014 Rim (L-zone) mean log-expression, '
             'z across arms', fontsize=7, pad=4)
for i in range(z.shape[0]):
    for j in range(z.shape[1]):
        ax.text(j, i, f'{pm.iloc[i, j]:.2f}', ha='center', va='center', fontsize=5,
                color='white' if abs(z.iloc[i, j]) > 1.3 else 'black')
cb = fig.colorbar(im, cax=fig.add_axes([0.968, 0.715, 0.010, 0.17]))
cb.set_label('z-score', fontsize=6); cb.ax.tick_params(labelsize=5.5)

# --- 4C rebuilt -----------------------------------------------------------
axC = [fig.add_subplot(gs[1:3, k]) for k in range(5)]
axC[0].text(-0.42, 1.05, 'C', transform=axC[0].transAxes,
            fontsize=10, fontweight='bold', va='top')
draw_C(axC, distribution=False)

fig.text(0.005, 0.985,
         'Figure 4.  Differential gene expression and pathway enrichment '
         '(PTPRC+ Rim immune spots)', fontsize=8, fontweight='bold', va='top')
footnote(fig, 0.062, False)
save(fig, 'Figure4_DEG_pathway_REGENERATED')

# ---------------------------------------------------------------- table ----
rows = []
for col, name, _ in PATHWAYS:
    for g in present:
        s = stat[col].loc[g]
        rows.append(dict(pathway_id=col, pathway=name.replace('\n', ' ').strip(),
                         arm=GROUP_LABEL[g],
                         n_spots=int(s['n']), mean=s['mean'], sem=s['sem'],
                         median=s['median']))
src = pd.DataFrame(rows)
src.to_csv(OUT / 'Figure4C_source_data.csv', index=False,
           encoding='utf-8-sig', float_format='%.6f')
print('  [saved] Figure4C_source_data.csv')
print('\n[done]', OUT)
