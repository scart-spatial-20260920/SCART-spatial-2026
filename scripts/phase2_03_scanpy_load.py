# ============================================================================
# phase2_03_scanpy_load.py
# Load the exported .h5ad from Phase 2 Step 2, verify metadata integrity,
# regenerate Fig 1 colored by actual cell type (not zone), and write
# per-group / per-zone summary tables.
#
# Prereq: the GEO h5ad is available under $SCART_DATA.
# ============================================================================

import os
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=120, facecolor="white")

DATA    = Path(os.environ.get("SCART_DATA", "./data"))
RESULTS = Path(os.environ.get("RESULTS", "./results"))
H5AD = DATA / "SCART_spatial_annotated.h5ad"   # name of the GEO-deposited AnnData
if not H5AD.exists():
    H5AD = DATA / "data_labeled.h5ad"          # original working-file name
FIG  = RESULTS / "figures"
QC   = RESULTS / "qc"
FIG.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

if not H5AD.exists():
    print(f"ERROR: {H5AD} not found. Download the GEO data and set SCART_DATA to that directory.")
    sys.exit(1)

print(f"[load] {H5AD}")
adata = sc.read_h5ad(H5AD)
print(adata)

print("\n[obs columns]"); print(adata.obs.columns.tolist())
print("\n[obs head]")
print(adata.obs.head(3).to_string())

# identify likely cell type column — prefer annotation columns over orig.ident
preferred = ["CellAnnotation", "cell_type", "celltype", "singleR_Ai", "seurat_clusters"]
ct_candidates = [c for c in preferred if c in adata.obs.columns]
if not ct_candidates:
    ct_candidates = [c for c in adata.obs.columns
                     if any(k in c.lower() for k in ["celltype","cell_type","cluster","annotation","singler"])]
print(f"\n[qc] cell type candidates: {ct_candidates}")

# identify spatial basis
if "spatial" not in adata.obsm and "X_spatial" not in adata.obsm:
    if {"spatial_x","spatial_y"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[["spatial_x","spatial_y"]].to_numpy(dtype=float)
        print("[qc] built adata.obsm['spatial'] from obs.spatial_x / spatial_y")

# phase2 labels present?
for c in ("phase2_group", "phase2_zone"):
    if c in adata.obs.columns:
        print(f"\n[qc] {c} counts:")
        print(adata.obs[c].value_counts(dropna=False))

# summary cross-tabulation: group x zone x cell type (if cell type exists)
if ct_candidates:
    ct_col = ct_candidates[0]
    print(f"\n[summary] using cell type column: {ct_col}")
    tab = adata.obs.groupby(["phase2_group","phase2_zone", ct_col]).size().unstack(fill_value=0)
    tab.to_csv(QC / "phase2_celltype_by_groupzone.csv")
    print(f"[write] {QC / 'phase2_celltype_by_groupzone.csv'}")

# spatial plot per group colored by cell type
if "spatial" in adata.obsm and ct_candidates:
    ct_col = ct_candidates[0]
    groups = ["control","SBRT","SCART1","SCART3"]
    fig, axes = plt.subplots(1, 4, figsize=(24, 6))
    for ax, g in zip(axes, groups):
        sub = adata[adata.obs["phase2_group"] == g].copy()
        if sub.n_obs == 0:
            ax.set_title(f"{g} (empty)"); continue
        x = sub.obsm["spatial"][:, 0]
        y = sub.obsm["spatial"][:, 1]
        labels = sub.obs[ct_col].astype(str).values
        uniq = sorted(pd.unique(labels))
        cmap = plt.get_cmap("tab20", max(len(uniq), 3))
        color_map = {l: cmap(i) for i, l in enumerate(uniq)}
        cols = [color_map[l] for l in labels]
        ax.scatter(x, y, s=1.2, c=cols, rasterized=True, alpha=0.85)
        ax.set_title(f"{g}  n={sub.n_obs:,}")
        ax.set_aspect("equal", adjustable="datalim")
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
    # legend from the full object
    all_labels = sorted(pd.unique(adata.obs[ct_col].astype(str)))
    cmap = plt.get_cmap("tab20", max(len(all_labels), 3))
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=cmap(i), label=l) for i, l in enumerate(all_labels)]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(all_labels), 8),
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.05))
    fig.suptitle("SCART Phase 2 — Fig 1C: spatial cell-type map by group", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG / "phase2_fig1C_spatial_celltype.png", dpi=300, bbox_inches="tight")
    plt.savefig(FIG / "phase2_fig1C_spatial_celltype.pdf", bbox_inches="tight")
    plt.close()
    print(f"[fig] wrote phase2_fig1C_spatial_celltype.png/.pdf")

print("[done] scanpy load and basic summary complete.")
