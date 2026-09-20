# ============================================================================
# phase2_08_morans_i.py
# Moran's I spatial autocorrelation per group on the 20-gene immune panel.
#
# Hypothesis: SCART3 has more spatially organized immune gene programs than SBRT.
# If true, SCART3 will show higher Moran's I for immune-related genes.
#
# Outputs:
#   phase2_qc/morans_i_immune_panel.csv
#   phase2_figures/phase2_figS_morans_i_panel.png
# ============================================================================

import os
from pathlib import Path
import sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

print = lambda *a, **k: (sys.__stdout__.write(" ".join(str(x) for x in a) + "\n"),
                          sys.__stdout__.flush())[0]

DATA    = Path(os.environ.get("SCART_DATA", "./data"))
RESULTS = Path(os.environ.get("RESULTS", "./results"))
H5AD = DATA / "SCART_spatial_annotated.h5ad"   # name of the GEO-deposited AnnData
if not H5AD.exists():
    H5AD = DATA / "data_labeled.h5ad"          # original working-file name
FIG  = RESULTS / "figures"
FIG.mkdir(parents=True, exist_ok=True)
QC   = RESULTS / "qc"
QC.mkdir(parents=True, exist_ok=True)

print("[load]")
adata = sc.read_h5ad(H5AD)
print(f"  {adata.n_obs:,} × {adata.n_vars:,}")

GROUPS = ["control", "SBRT", "SCART1", "SCART3"]

# Phase 1 immune panel
PANEL = ["H2-Ab1","H2-Aa","Ccr7","Cd86","Cd74",
         "Isg15","Ifit1","Ifit3","Cxcl10","Cxcl9","Ifng",
         "Foxp3","Il10","Tgfb1","Pdcd1","Lag3",
         "Nos2","Il1b","Arg1","Mrc1"]
PANEL_PRESENT = [g for g in PANEL if g in adata.var_names]
print(f"  panel: {len(PANEL_PRESENT)}/{len(PANEL)} genes present")
PANEL_MISSING = set(PANEL) - set(PANEL_PRESENT)
if PANEL_MISSING:
    print(f"  missing: {PANEL_MISSING}")

# Add a few extra genes of interest for spatial signal
EXTRA = ["C1qa","C1qb","C1qc","B2m","Apoe","Lyz2","Ctss","S100a4",
         "Ccl2","Ccl5","Cd3e","Cd4","Cd8a","Gzmb","Prf1","Nkg7"]
EXTRA_PRESENT = [g for g in EXTRA if g in adata.var_names]
ALL_GENES = list(set(PANEL_PRESENT + EXTRA_PRESENT))
print(f"  total genes for Moran's I: {len(ALL_GENES)}")

# ---------- compute Moran's I per group ----------
results = {}
for g in GROUPS:
    print(f"\n[moransI] {g} ...")
    sub = adata[adata.obs["phase2_group"] == g].copy()
    if sub.n_obs < 50:
        print(f"  too few cells ({sub.n_obs}), skip")
        continue
    # subset to genes of interest (faster)
    sub_genes = sub[:, [g for g in ALL_GENES if g in sub.var_names]].copy()
    # spatial neighbors
    sq.gr.spatial_neighbors(sub_genes, coord_type="generic", n_neighs=10)
    # Moran's I
    sq.gr.spatial_autocorr(sub_genes, mode="moran", genes=sub_genes.var_names.tolist(), n_jobs=1)
    mi = sub_genes.uns["moranI"]
    mi["group"] = g
    results[g] = mi
    print(f"  n_cells={sub.n_obs:,}  genes={mi.shape[0]}  median I={mi['I'].median():.4f}")
    # top 5
    print("  top 5 by Moran's I:")
    print(mi.nlargest(5, "I")[["I","pval_norm"]].to_string())

# combine
all_mi = pd.concat(results.values(), ignore_index=False)
all_mi.to_csv(QC / "morans_i_immune_panel.csv")
print(f"\n[write] morans_i_immune_panel.csv ({len(all_mi)} rows)")

# ---------- heatmap: Moran's I per gene × group ----------
pivot = all_mi.reset_index().pivot_table(index="index", columns="group", values="I")
pivot = pivot.reindex(columns=GROUPS)

# sort by SCART3 Moran's I descending
if "SCART3" in pivot.columns:
    pivot = pivot.sort_values("SCART3", ascending=False)

fig, ax = plt.subplots(figsize=(6, max(8, len(pivot)*0.22)))
im = ax.imshow(pivot.values, cmap="YlOrRd", aspect="auto", vmin=0)
ax.set_xticks(range(pivot.shape[1]))
ax.set_xticklabels(pivot.columns, fontsize=10)
ax.set_yticks(range(pivot.shape[0]))
ax.set_yticklabels(pivot.index, fontsize=8)
plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="Moran's I")
ax.set_title("Moran's I spatial autocorrelation (immune panel)\nhigher = more spatially structured")
plt.tight_layout()
plt.savefig(FIG / "phase2_figS_morans_i_panel.png", dpi=250, bbox_inches="tight")
plt.savefig(FIG / "phase2_figS_morans_i_panel.pdf", bbox_inches="tight")
plt.close()
print(f"[fig] phase2_figS_morans_i_panel.png/.pdf")

# ---------- summary: mean Moran's I per group ----------
summary = pivot.mean(axis=0)
print("\n[summary] mean Moran's I per group:")
for g in GROUPS:
    if g in summary.index:
        print(f"  {g}: {summary[g]:.4f}")
best = summary.idxmax()
print(f"\n  HIGHEST mean spatial structure: {best} ({summary[best]:.4f})")

print("\n[done]")
