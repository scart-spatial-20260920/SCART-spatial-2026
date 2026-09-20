"""
phase2_10_journal_figures.py
Journal-submission figure set for spatial transcriptomics (SCART vs SBRT, mouse HCC).

Specs locked to manuscript requirements:
  - Vector PDF + 300-DPI PNG + 600-DPI TIFF for each figure
  - Arial 6-8 pt (no overlap, no blur on Fig 3A heatmap)
  - Zone palette  STV (H) = red, TTV (M) = orange, Rim (L) = blue
  - 14-type immune palette: tab20 (Nature-style)
  - 5 Hallmark pathways: IFN-gamma, Inflammatory, TNFA/NFkB, IL6/JAK/STAT3, Complement

Output:   phase2_figures/journal/
"""

import os
from pathlib import Path
import sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.colors import TwoSlopeNorm
import scanpy as sc

# --- paths ---------------------------------------------------------------
DATA    = Path(os.environ.get("SCART_DATA", "./data"))
RESULTS = Path(os.environ.get("RESULTS", "./results"))
H5AD = DATA / "SCART_spatial_annotated.h5ad"   # name of the GEO-deposited AnnData
if not H5AD.exists():
    H5AD = DATA / "data_labeled.h5ad"          # original working-file name
QC   = RESULTS / "qc"
QC.mkdir(parents=True, exist_ok=True)
OUT  = RESULTS / "figures" / "journal"
OUT.mkdir(parents=True, exist_ok=True)

# --- style ---------------------------------------------------------------
plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "legend.title_fontsize": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "lines.linewidth": 0.8,
    "pdf.fonttype": 42,   # editable text in PDF
    "ps.fonttype": 42,
})

GROUPS = ["control", "SBRT", "SCART1", "SCART3"]
GROUP_COLOR = {"control": "#7f7f7f", "SBRT": "#4575b4", "SCART1": "#fdae61", "SCART3": "#d73027"}
ZONE_COLOR  = {"H": "#d62728", "M": "#ff7f0e", "L": "#1f77b4", "X": "#e8e8e8"}
ZONE_LABEL  = {"H": "STV (core)", "M": "TTV (junction)", "L": "Rim (periphery)", "X": "Outside zones"}

def save(fig, name):
    """PNG + PDF + TIFF (600 dpi)."""
    fig.savefig(OUT / f"{name}.png", dpi=300, facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", facecolor="white")
    fig.savefig(OUT / f"{name}.tiff", dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"[saved] {name}.{{png,pdf,tiff}}")

def panel_label(ax, letter, x=-0.10, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")

# --- load ---------------------------------------------------------------
print("[load]", H5AD)
adata = sc.read_h5ad(H5AD)

# merge tissue-class spelling variants
ai = adata.obs["singleR_Ai"].astype(str).replace({
    "Endothelial_cell": "Endothelial",
    "Endothelial Cell": "Endothelial",
})
adata.obs["tissue_class"] = pd.Categorical(ai)

fine_col = "GPT_PTPRC_and_PTPRCAP"
imm_mask = adata.obs[fine_col].astype(str) != "Undefined"
immune_levels = sorted(adata.obs.loc[imm_mask, fine_col].unique().tolist())
tab20 = plt.get_cmap("tab20", max(len(immune_levels), 3))
imm_cmap = {lv: tab20(i) for i, lv in enumerate(immune_levels)}

print(f"  cells={adata.n_obs:,}  immune={imm_mask.sum():,}  types={len(immune_levels)}")

# ===========================================================================
# FIGURE 1.  Spatial overview
#   A-D: 4 groups × all cells × STV/TTV/Rim zones (red/orange/blue)
#   E  : PTPRC+ immune cells (14 types) on tissue coordinates, 4 groups
# ===========================================================================
print("\n[Figure 1] spatial overview")
fig = plt.figure(figsize=(7.2, 4.2))   # 2-column Nature width
gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.32, wspace=0.10,
                       left=0.04, right=0.84, top=0.86, bottom=0.05)

for i, g in enumerate(GROUPS):
    # --- row 1: zone maps ---
    ax = fig.add_subplot(gs[0, i])
    sub = adata[adata.obs["phase2_group"] == g]
    xy  = sub.obsm["spatial"]
    z   = sub.obs["phase2_zone"].astype(str).values
    for zn in ["X", "L", "M", "H"]:   # background first, then L/M/H above
        m = z == zn
        if m.sum() == 0: continue
        ax.scatter(xy[m, 0], xy[m, 1], s=0.25,
                   c=ZONE_COLOR[zn], alpha=0.8, rasterized=True, linewidths=0)
    ax.set_aspect("equal", adjustable="datalim"); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(f"{g}\nn = {sub.n_obs:,}", fontsize=7, pad=2)
    if i == 0: panel_label(ax, "A")
    elif i == 1: panel_label(ax, "B")
    elif i == 2: panel_label(ax, "C")
    elif i == 3: panel_label(ax, "D")

    # --- row 2: immune 14-type maps ---
    ax2 = fig.add_subplot(gs[1, i])
    sub_all = adata[adata.obs["phase2_group"] == g]
    ax2.scatter(sub_all.obsm["spatial"][:, 0], sub_all.obsm["spatial"][:, 1],
                s=0.15, c="#dddddd", alpha=0.5, rasterized=True, linewidths=0)
    imm_sub = adata[(adata.obs["phase2_group"] == g) & imm_mask]
    if imm_sub.n_obs > 0:
        cols = np.array([imm_cmap[v] for v in imm_sub.obs[fine_col].astype(str)])
        ax2.scatter(imm_sub.obsm["spatial"][:, 0], imm_sub.obsm["spatial"][:, 1],
                    s=2.5, c=cols, alpha=0.95, rasterized=True,
                    edgecolors="black", linewidths=0.08)
    ax2.set_aspect("equal", adjustable="datalim"); ax2.invert_yaxis()
    ax2.set_xticks([]); ax2.set_yticks([])
    for sp in ax2.spines.values(): sp.set_visible(False)
    ax2.set_title(f"{g}\nPTPRC+  n = {imm_sub.n_obs}", fontsize=7, pad=2)
    if i == 0: panel_label(ax2, "E")

# legends on the right edge
zone_handles = [Patch(fc=ZONE_COLOR[k], ec="none", label=ZONE_LABEL[k])
                for k in ["H", "M", "L", "X"]]
fig.legend(handles=zone_handles, loc="upper left",
           bbox_to_anchor=(0.85, 0.92), frameon=False,
           title="Spatial zone", fontsize=6, title_fontsize=7)

imm_handles = [Patch(fc=imm_cmap[lv], ec="none", label=lv) for lv in immune_levels]
fig.legend(handles=imm_handles, loc="upper left",
           bbox_to_anchor=(0.85, 0.48), frameon=False,
           title="Immune (PTPRC+) subtype", fontsize=5.5, title_fontsize=7,
           ncol=1, handlelength=1, handletextpad=0.4, labelspacing=0.3)

fig.text(0.02, 0.97, "Figure 1.  Spatial overview", fontsize=8,
         fontweight="bold")
save(fig, "Figure1_spatial_overview")

# ===========================================================================
# FIGURE 2.  Tissue architecture by broad cell class (singleR)
# ===========================================================================
print("\n[Figure 2] tissue architecture")
broad_order = [
    "Endothelial", "Fibroblasts", "Hepatocyte",
    "Immune_Cell", "Myeloid_cell", "Macrophage_Kupffer",
    "Proliferating Cells", "Unclassified",
]
broad_present = [b for b in broad_order
                 if b in adata.obs["tissue_class"].astype(str).unique()]
tab10 = plt.get_cmap("tab10", max(len(broad_present), 3))
broad_cmap = {lv: tab10(i) for i, lv in enumerate(broad_present)}

fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.6))
for i, (ax, g) in enumerate(zip(axes, GROUPS)):
    sub = adata[adata.obs["phase2_group"] == g]
    xy  = sub.obsm["spatial"]
    tc  = sub.obs["tissue_class"].astype(str).values
    cols = np.array([broad_cmap.get(v, "#cccccc") for v in tc])
    ax.scatter(xy[:, 0], xy[:, 1], s=0.4, c=cols, alpha=0.85,
               rasterized=True, linewidths=0)
    ax.set_aspect("equal", adjustable="datalim"); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(f"{g}\nn = {sub.n_obs:,}", fontsize=7, pad=2)
    if i == 0: panel_label(ax, "A", x=-0.05)

handles = [Patch(fc=broad_cmap[lv], ec="none", label=lv) for lv in broad_present]
fig.legend(handles=handles, loc="lower center",
           bbox_to_anchor=(0.5, -0.07), ncol=4, frameon=False, fontsize=6)
fig.text(0.02, 0.98, "Figure 2.  Tissue architecture by broad cell class",
         fontsize=8, fontweight="bold")
plt.subplots_adjust(left=0.04, right=0.98, top=0.78, bottom=0.12, wspace=0.10)
save(fig, "Figure2_tissue_architecture")

# ===========================================================================
# FIGURE 3.  Spatial immune ecology
#   A: neighborhood enrichment z-scores (broad classes, 4 groups, KNN=10, 200 perm)
#   B: SBRT (left) vs SCART3 (right) rim immune cell maps
# ===========================================================================
print("\n[Figure 3] spatial immune ecology")
z_all   = pd.read_csv(QC / "nhood_enrichment_zscores.csv")
z_broad = z_all[z_all["annotation"] == "broad"].copy()
cats = sorted(set(z_broad["cluster_a"]) | set(z_broad["cluster_b"]))

fig = plt.figure(figsize=(7.2, 6.6))
gs = gridspec.GridSpec(2, 4, figure=fig,
                       height_ratios=[1.0, 1.0],
                       hspace=0.42, wspace=0.10,
                       left=0.10, right=0.95, top=0.94, bottom=0.06)

# 3A — 4 heatmaps + shared colorbar
mats = {}
for g in GROUPS:
    sub = z_broad[z_broad["group"] == g]
    if sub.empty: continue
    piv = sub.pivot_table(index="cluster_a", columns="cluster_b",
                          values="zscore", aggfunc="mean")
    mats[g] = piv.reindex(index=cats, columns=cats)

vmax = max(np.nanmax(np.abs(m.values)) for m in mats.values())
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

for j, g in enumerate(GROUPS):
    ax = fig.add_subplot(gs[0, j])
    if j == 0: panel_label(ax, "A", x=-0.45)
    m = mats[g].fillna(0)
    im = ax.imshow(m.values, cmap="RdBu_r", norm=norm, aspect="equal")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, rotation=55, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(cats)))
    ax.set_yticklabels(cats if j == 0 else [], fontsize=5.5)
    ax.set_title(g, fontsize=8, pad=4)
    ax.tick_params(axis="both", which="major", pad=1, length=2)
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_linewidth(0.4)

# colorbar
cbar_ax = fig.add_axes([0.96, 0.55, 0.012, 0.30])
cb = fig.colorbar(im, cax=cbar_ax)
cb.set_label("Neighborhood enrichment z-score", fontsize=6)
cb.ax.tick_params(labelsize=5.5)

# 3B — SBRT vs SCART3 rim immune detail
for j, g in enumerate(["SBRT", "SCART3"]):
    ax = fig.add_subplot(gs[1, j*2:(j+1)*2])
    if j == 0: panel_label(ax, "B", x=-0.05)
    grp = adata.obs["phase2_group"] == g
    bg = adata[grp]
    ax.scatter(bg.obsm["spatial"][:, 0], bg.obsm["spatial"][:, 1],
               s=0.2, c="#eaeaea", alpha=0.55, rasterized=True, linewidths=0)
    l_imm = adata[grp & (adata.obs["phase2_zone"] == "L") & imm_mask]
    if l_imm.n_obs > 0:
        cols = np.array([imm_cmap[v] for v in l_imm.obs[fine_col].astype(str)])
        ax.scatter(l_imm.obsm["spatial"][:, 0], l_imm.obsm["spatial"][:, 1],
                   s=8, c=cols, alpha=0.95, rasterized=True,
                   edgecolors="black", linewidths=0.15)
    ax.set_aspect("equal", adjustable="datalim"); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(f"{g} — Rim immune  n = {l_imm.n_obs}", fontsize=8, pad=2)

# immune legend below row 2 (compact 2-column)
imm_legend_ax = fig.add_axes([0.05, -0.02, 0.9, 0.03])
imm_legend_ax.axis("off")
imm_legend_ax.legend(handles=imm_handles, loc="center", ncol=7, frameon=False,
                      fontsize=5.5, handlelength=1, handletextpad=0.3,
                      columnspacing=0.8)

fig.text(0.02, 0.97, "Figure 3.  Spatial immune ecology", fontsize=8,
         fontweight="bold")
fig.text(0.02, 0.94,
         "Neighborhood enrichment (Squidpy, KNN=10, 200 permutations); "
         "diagonal = self-aggregation.",
         fontsize=6, style="italic")
save(fig, "Figure3_spatial_immune_ecology")

# ===========================================================================
# FIGURE 4.  DEG and pathway enrichment
#   A: volcano SCART3_L vs SBRT_L (PTPRC+, n=393 vs 230)
#   B: phase-1 immune gene panel (z-score across groups, L-zone)
#   C: 5 Hallmark pathways (IFN-g, inflam, TNFA, IL6/JAK/STAT3, complement)
# ===========================================================================
print("\n[Figure 4] DEG + pathways")
deg = pd.read_csv(QC / "deg_SCART3L_vs_SBRTL.csv")
panel_means = pd.read_csv(QC / "phase1_panel_means_in_L_zone.csv", index_col=0)
pw_df = pd.read_csv(QC / "pathway_aucell_per_cell.csv", index_col=0)

# count the n used by the DEG run
n_scart3_rim = int(((adata.obs["phase2_group"]=="SCART3") &
                    (adata.obs["phase2_zone"]=="L") & imm_mask).sum())
n_sbrt_rim   = int(((adata.obs["phase2_group"]=="SBRT") &
                    (adata.obs["phase2_zone"]=="L") & imm_mask).sum())

fig = plt.figure(figsize=(7.2, 7.2))
gs  = gridspec.GridSpec(3, 5, figure=fig,
                        height_ratios=[1.2, 1.2, 1.0],
                        hspace=0.55, wspace=0.55,
                        left=0.09, right=0.97, top=0.94, bottom=0.07)

# --- 4A volcano ---
ax = fig.add_subplot(gs[0, 0:2])
panel_label(ax, "A", x=-0.18)
v = deg.copy()
v["mlp"] = -np.log10(v["pvals_adj"].clip(lower=1e-300))
v.loc[v["mlp"] > 50, "mlp"] = 50
sig = (v["pvals_adj"] < 0.05) & (v["logfoldchanges"].abs() > 0.5)
ax.scatter(v.loc[~sig, "logfoldchanges"], v.loc[~sig, "mlp"],
           s=2.5, c="#cccccc", alpha=0.45, rasterized=True, linewidths=0)
ax.scatter(v.loc[sig & (v["logfoldchanges"]>0), "logfoldchanges"],
           v.loc[sig & (v["logfoldchanges"]>0), "mlp"],
           s=5, c="#d73027", alpha=0.85,
           label=f"↑ SCART3 (n={n_scart3_rim})", rasterized=True, linewidths=0)
ax.scatter(v.loc[sig & (v["logfoldchanges"]<0), "logfoldchanges"],
           v.loc[sig & (v["logfoldchanges"]<0), "mlp"],
           s=5, c="#4575b4", alpha=0.85,
           label=f"↑ SBRT (n={n_sbrt_rim})", rasterized=True, linewidths=0)
top_up   = v.loc[sig & (v["logfoldchanges"]>0)].nlargest(8, "logfoldchanges")
top_down = v.loc[sig & (v["logfoldchanges"]<0)].nsmallest(8, "logfoldchanges")
for _, r in pd.concat([top_up, top_down]).iterrows():
    ax.annotate(r["names"], (r["logfoldchanges"], r["mlp"]),
                fontsize=5, alpha=0.85,
                xytext=(2, 1), textcoords="offset points")
ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.4)
ax.axvline(0.5,  color="grey", ls=":", lw=0.4)
ax.axvline(-0.5, color="grey", ls=":", lw=0.4)
ax.set_xlabel("log₂ fold-change (SCART3 / SBRT)", fontsize=7)
ax.set_ylabel("−log₁₀ (adjusted p)", fontsize=7)
ax.set_title("DEG  —  SCART3 Rim vs SBRT Rim  (PTPRC+ immune)",
             fontsize=7.5, pad=3)
ax.legend(fontsize=6, loc="upper left", frameon=False)

# --- 4B immune panel heatmap ---
ax = fig.add_subplot(gs[0, 2:5])
panel_label(ax, "B", x=-0.10)
pm = panel_means.copy()
pm.columns = [c.replace("_L", "") for c in pm.columns]
pm = pm.reindex(columns=GROUPS)
z = pm.sub(pm.mean(axis=1), axis=0).div(pm.std(axis=1).replace(0, 1), axis=0)
im = ax.imshow(z.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
ax.set_xticks(range(z.shape[1]))
ax.set_xticklabels(z.columns, fontsize=7)
ax.set_yticks(range(z.shape[0]))
ax.set_yticklabels(z.index, fontsize=6)
ax.tick_params(axis="y", which="major", pad=1, length=2)
ax.set_title("Phase-1 immune gene panel  —  Rim (L-zone)  mean log-expression  /  z across groups",
             fontsize=7.5, pad=3)
for i in range(z.shape[0]):
    for j in range(z.shape[1]):
        val = pm.iloc[i, j]
        ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5,
                color="white" if abs(z.iloc[i, j]) > 1.3 else "black")
cbax = fig.add_axes([0.965, 0.62, 0.010, 0.20])
cb = fig.colorbar(im, cax=cbax); cb.set_label("z-score", fontsize=6)
cb.ax.tick_params(labelsize=5.5)

# --- 4C pathway bars (5 pathways) ---
target = [
    ("HALLMARK_INTERFERON_GAMMA_RESPONSE", "IFN-γ response"),
    ("HALLMARK_INFLAMMATORY_RESPONSE",     "Inflammatory response"),
    ("HALLMARK_TNFA_SIGNALING_VIA_NFKB",   "TNFα / NF-κB"),
    ("HALLMARK_IL6_JAK_STAT3_SIGNALING",   "IL6 / JAK / STAT3"),
    ("HALLMARK_COMPLEMENT",                "Complement"),
]
agg_mean = pw_df.groupby("phase2_group").mean(numeric_only=True)
agg_sem  = pw_df.groupby("phase2_group").sem(numeric_only=True)

for k, (col, lab) in enumerate(target):
    ax = fig.add_subplot(gs[1:3, k])
    if k == 0: panel_label(ax, "C", x=-0.30)
    if col not in agg_mean.columns:
        ax.axis("off"); continue
    vals = agg_mean.loc[[g for g in GROUPS if g in agg_mean.index], col]
    errs = agg_sem.loc[vals.index,  col]
    ax.bar(range(len(vals)), vals.values, yerr=errs.values,
           color=[GROUP_COLOR[g] for g in vals.index],
           edgecolor="black", linewidth=0.5, capsize=2,
           error_kw=dict(elinewidth=0.5, capthick=0.5))
    ax.axhline(0, color="black", lw=0.4)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(vals.index, rotation=40, ha="right", fontsize=6)
    ax.set_title(lab, fontsize=7, pad=3)
    if k == 0: ax.set_ylabel("Pathway score  (mean ± SEM)", fontsize=6.5)

fig.text(0.02, 0.97,
         "Figure 4.  Differential gene expression and pathway enrichment  (PTPRC+ Rim immune cells)",
         fontsize=8, fontweight="bold")
save(fig, "Figure4_DEG_pathway")

# ===========================================================================
# manifest
# ===========================================================================
print("\n[done]  journal figures in:")
print(f"  {OUT}")
for f in sorted(OUT.glob("Figure*.*")):
    sz = f.stat().st_size / 1e6
    print(f"   {f.name:48s}  {sz:6.2f} MB")
