# ============================================================================
# phase2_04_spatial_celltype.py
# Proper spatial cell-type maps for the manuscript:
#   - Fig 1C(i): tissue context — all cells colored by broad class (singleR_Ai)
#   - Fig 1C(ii): immune subset colored by 14-type fine annotation (GPT_PTPRC_and_PTPRCAP)
#   - Zoom insets for the Rim (L-zone) of SBRT and SCART3 — the key story
# ============================================================================

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import scanpy as sc

sc.settings.verbosity = 1

DATA    = Path(os.environ.get("SCART_DATA", "./data"))
RESULTS = Path(os.environ.get("RESULTS", "./results"))
H5AD = DATA / "SCART_spatial_annotated.h5ad"   # name of the GEO-deposited AnnData
if not H5AD.exists():
    H5AD = DATA / "data_labeled.h5ad"          # original working-file name
FIG  = RESULTS / "figures"
FIG.mkdir(parents=True, exist_ok=True)
QC   = RESULTS / "qc"
QC.mkdir(parents=True, exist_ok=True)

print(f"[load] {H5AD}")
adata = sc.read_h5ad(H5AD)
print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

GROUPS = ["control", "SBRT", "SCART1", "SCART3"]
ZONE_COLOR = {"H": "#d73027", "M": "#fc8d59", "L": "#4575b4", "X": "#dddddd"}

# -------- broad tissue class (singleR_Ai) --------
# normalize some duplicate labels
ai = adata.obs["singleR_Ai"].astype(str).copy()
ai = ai.replace({"Endothelial_cell": "Endothelial",
                 "Endothelial Cell": "Endothelial",
                 "Endothelial": "Endothelial"})
adata.obs["tissue_class"] = ai.astype("category")
broad_levels = adata.obs["tissue_class"].value_counts().index.tolist()
print(f"[broad] tissue_class levels: {broad_levels}")

# -------- 14-type fine immune annotation --------
fine_col = "GPT_PTPRC_and_PTPRCAP"
fine = adata.obs[fine_col].astype(str)
immune_mask = fine != "Undefined"
print(f"[fine] immune cells (PTPRC+): {immune_mask.sum():,} / {adata.n_obs:,}")

# count per zone per group using the ground-truth 12-region label
tab = (adata.obs.loc[immune_mask]
       .groupby([fine_col, "P_group_all_HML"], observed=True)
       .size().unstack(fill_value=0))
tab.to_csv(QC / "immune_celltype_by_region.csv")
print("[write] immune_celltype_by_region.csv")

# =============================================================================
# Figure 1C panel (i): all cells by tissue class, 4-panel
# =============================================================================
fig, axes = plt.subplots(1, 4, figsize=(24, 6))
tab20 = plt.get_cmap("tab10", max(len(broad_levels), 3))
broad_color_map = {lv: tab20(i) for i, lv in enumerate(broad_levels)}

for ax, g in zip(axes, GROUPS):
    sub = adata[adata.obs["phase2_group"] == g]
    if sub.n_obs == 0:
        ax.set_title(f"{g} (empty)"); continue
    xy = sub.obsm["spatial"]
    cols = np.array([broad_color_map[v] for v in sub.obs["tissue_class"].astype(str)])
    ax.scatter(xy[:,0], xy[:,1], s=0.8, c=cols, rasterized=True, alpha=0.85)
    ax.set_title(f"{g}  n={sub.n_obs:,}", fontsize=12)
    ax.set_aspect("equal", adjustable="datalim")
    ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])

handles = [Patch(facecolor=broad_color_map[lv], label=f"{lv}") for lv in broad_levels]
fig.legend(handles=handles, loc="lower center", ncol=min(len(broad_levels), 6),
           frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.04))
fig.suptitle("SCART Phase 2 — Fig 1C(i): tissue class map (singleR broad annotation)", fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(FIG / "phase2_fig1C_i_tissue_class.png", dpi=300, bbox_inches="tight")
plt.savefig(FIG / "phase2_fig1C_i_tissue_class.pdf", bbox_inches="tight")
plt.close()
print("[fig] phase2_fig1C_i_tissue_class.png/.pdf")

# =============================================================================
# Figure 1C panel (ii): immune subset 14-type, 4-panel with tissue background
# =============================================================================
fig, axes = plt.subplots(1, 4, figsize=(24, 6))

immune_levels = sorted(adata.obs.loc[immune_mask, fine_col].unique().tolist())
tab20 = plt.get_cmap("tab20", max(len(immune_levels), 3))
immune_color_map = {lv: tab20(i) for i, lv in enumerate(immune_levels)}

for ax, g in zip(axes, GROUPS):
    grp_mask = adata.obs["phase2_group"] == g
    sub_all = adata[grp_mask]
    # tissue background
    bg_xy = sub_all.obsm["spatial"]
    ax.scatter(bg_xy[:,0], bg_xy[:,1], s=0.3, c="#eeeeee", rasterized=True, alpha=0.5)
    # immune cells on top
    imm_sub = adata[grp_mask & immune_mask]
    if imm_sub.n_obs > 0:
        xy = imm_sub.obsm["spatial"]
        cols = np.array([immune_color_map[v] for v in imm_sub.obs[fine_col].astype(str)])
        ax.scatter(xy[:,0], xy[:,1], s=8, c=cols, rasterized=True, alpha=0.9,
                   edgecolors="black", linewidths=0.15)
    ax.set_title(f"{g}   immune n={imm_sub.n_obs:,}", fontsize=12)
    ax.set_aspect("equal", adjustable="datalim")
    ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])

handles = [Patch(facecolor=immune_color_map[lv], label=lv) for lv in immune_levels]
fig.legend(handles=handles, loc="lower center", ncol=5,
           frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.08))
fig.suptitle("SCART Phase 2 — Fig 1C(ii): PTPRC+ immune cells colored by 14-type annotation", fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(FIG / "phase2_fig1C_ii_immune_14type.png", dpi=300, bbox_inches="tight")
plt.savefig(FIG / "phase2_fig1C_ii_immune_14type.pdf", bbox_inches="tight")
plt.close()
print("[fig] phase2_fig1C_ii_immune_14type.png/.pdf")

# =============================================================================
# Figure 1C panel (iii): SBRT_L vs SCART3_L zoom — the rim showdown
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 7))
for ax, g in zip(axes, ["SBRT", "SCART3"]):
    grp_mask = (adata.obs["phase2_group"] == g) & (adata.obs["phase2_zone"] == "L")
    sub_all = adata[adata.obs["phase2_group"] == g]  # tissue context
    ax.scatter(sub_all.obsm["spatial"][:,0], sub_all.obsm["spatial"][:,1],
               s=0.4, c="#f2f2f2", rasterized=True)
    imm = adata[grp_mask & immune_mask]
    if imm.n_obs > 0:
        cols = np.array([immune_color_map[v] for v in imm.obs[fine_col].astype(str)])
        ax.scatter(imm.obsm["spatial"][:,0], imm.obsm["spatial"][:,1],
                   s=14, c=cols, alpha=0.95, edgecolors="black", linewidths=0.25, rasterized=True)
    # zoom to L-zone bounding box
    l_sub = adata[grp_mask]
    if l_sub.n_obs > 0:
        x = l_sub.obsm["spatial"][:,0]; y = l_sub.obsm["spatial"][:,1]
        ax.set_xlim(x.min()-200, x.max()+200)
        ax.set_ylim(y.max()+200, y.min()-200)  # inverted
    ax.set_title(f"{g}  Rim (L-zone) immune detail   n_immune={imm.n_obs:,}", fontsize=12)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([]); ax.set_yticks([])

fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=8,
           bbox_to_anchor=(0.5, -0.08))
fig.suptitle("SCART Phase 2 — Fig 3B: SBRT vs SCART3 rim immune ecology", fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(FIG / "phase2_fig3B_rim_sbrt_vs_scart3.png", dpi=300, bbox_inches="tight")
plt.savefig(FIG / "phase2_fig3B_rim_sbrt_vs_scart3.pdf", bbox_inches="tight")
plt.close()
print("[fig] phase2_fig3B_rim_sbrt_vs_scart3.png/.pdf")

print("[done]")
