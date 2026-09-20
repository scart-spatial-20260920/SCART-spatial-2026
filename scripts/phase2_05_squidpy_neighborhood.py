# ============================================================================
# phase2_05_squidpy_neighborhood.py
# Spatial neighborhood enrichment per group (Squidpy).
#
# Key question: does SCART3_L (rim) show stronger DC↔T-cell and inflam-mono↔TAM
# neighborhood enrichment than SBRT_L?
#
# Two runs:
#   (A) broad: tissue_class (all 58k cells) — denser graph, shows stromal context
#   (B) fine:  GPT_PTPRC_and_PTPRCAP (3,264 PTPRC+ cells) — the 14-type immune story
#
# Outputs:
#   phase2_figures/phase2_fig3A_nhood_broad_per_group.png
#   phase2_figures/phase2_fig3A_nhood_immune_per_group.png
#   phase2_figures/phase2_fig3B_scart3L_vs_sbrtL_immune_delta.png
#   phase2_qc/nhood_enrichment_zscores.csv
# ============================================================================

import os
from pathlib import Path
import sys
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

# unbuffered prints so we can watch progress in real time
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

print(f"[load] {H5AD}")
adata = sc.read_h5ad(H5AD)
print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# normalize tissue_class (duplicate Endothelial labels)
ai = adata.obs["singleR_Ai"].astype(str).copy()
ai = ai.replace({"Endothelial_cell": "Endothelial",
                 "Endothelial Cell": "Endothelial"})
adata.obs["tissue_class"] = ai.astype("category")

GROUPS = ["control", "SBRT", "SCART1", "SCART3"]

# ---------- helper: build spatial graph per-group and run nhood_enrichment ----------
def nhood_for_group(ad_all, group, cluster_key, radius=None, n_neighs=10):
    sub = ad_all[ad_all.obs["phase2_group"] == group].copy()
    if sub.n_obs < 30:
        return None, None
    # drop categories not present to avoid empty-cluster errors
    sub.obs[cluster_key] = sub.obs[cluster_key].astype(str).astype("category")
    sub.obs[cluster_key] = sub.obs[cluster_key].cat.remove_unused_categories()
    # spatial neighbors — KNN on spatial coords
    sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=n_neighs)
    # n_jobs=1 to avoid Windows multiprocessing deadlock
    sq.gr.nhood_enrichment(sub, cluster_key=cluster_key, seed=42, n_perms=200, n_jobs=1, show_progress_bar=False)
    zscore = sub.uns[f"{cluster_key}_nhood_enrichment"]["zscore"]
    count  = sub.uns[f"{cluster_key}_nhood_enrichment"]["count"]
    cats = sub.obs[cluster_key].cat.categories.tolist()
    z_df = pd.DataFrame(zscore, index=cats, columns=cats)
    c_df = pd.DataFrame(count,  index=cats, columns=cats)
    return z_df, c_df

# ============================================================================
# (A) Broad tissue class enrichment — all 58k cells
# ============================================================================
print("\n[A] Broad tissue class (singleR_Ai) — all cells")
broad_results = {}
for g in GROUPS:
    print(f"  running {g} ...")
    z, c = nhood_for_group(adata, g, "tissue_class", n_neighs=10)
    if z is not None:
        broad_results[g] = z
        print(f"    {z.shape[0]} clusters, max |z|={np.abs(z.values).max():.1f}")

# shared cluster ordering for plot
all_cats_broad = sorted(set().union(*[df.index for df in broad_results.values()]))
fig, axes = plt.subplots(1, 4, figsize=(22, 5.5))
vmax = max(np.abs(df.values).max() for df in broad_results.values())
for ax, g in zip(axes, GROUPS):
    df = broad_results.get(g)
    if df is None:
        ax.set_title(f"{g} (empty)"); ax.axis("off"); continue
    df = df.reindex(index=all_cats_broad, columns=all_cats_broad)
    im = ax.imshow(df.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(all_cats_broad)))
    ax.set_yticks(range(len(all_cats_broad)))
    ax.set_xticklabels(all_cats_broad, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(all_cats_broad, fontsize=8)
    ax.set_title(f"{g}   n={(adata.obs['phase2_group']==g).sum():,}", fontsize=11)
fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, label="nhood enrichment z-score")
fig.suptitle("SCART Phase 2 — Fig 3A(i): neighborhood enrichment by tissue class (broad)", fontsize=13, y=1.02)
plt.savefig(FIG / "phase2_fig3A_nhood_broad_per_group.png", dpi=250, bbox_inches="tight")
plt.savefig(FIG / "phase2_fig3A_nhood_broad_per_group.pdf", bbox_inches="tight")
plt.close()
print(f"  [fig] phase2_fig3A_nhood_broad_per_group.png/.pdf")

# ============================================================================
# (B) Fine immune 14-type enrichment — PTPRC+ cells only
# ============================================================================
print("\n[B] Fine immune 14-type (GPT_PTPRC_and_PTPRCAP) — PTPRC+ only")
fine_col = "GPT_PTPRC_and_PTPRCAP"
imm_mask = adata.obs[fine_col].astype(str) != "Undefined"
print(f"  {imm_mask.sum():,} immune cells")

imm = adata[imm_mask].copy()
fine_results = {}
for g in GROUPS:
    print(f"  running {g} ...")
    z, c = nhood_for_group(imm, g, fine_col, n_neighs=6)
    if z is not None:
        fine_results[g] = z
        print(f"    n_clusters={z.shape[0]}  max |z|={np.abs(z.values).max():.1f}")

all_cats_fine = sorted(set().union(*[df.index for df in fine_results.values()]))
fig, axes = plt.subplots(1, 4, figsize=(26, 7))
vmax = max(np.abs(df.values).max() for df in fine_results.values()) if fine_results else 5
for ax, g in zip(axes, GROUPS):
    df = fine_results.get(g)
    if df is None:
        ax.set_title(f"{g} (empty)"); ax.axis("off"); continue
    df = df.reindex(index=all_cats_fine, columns=all_cats_fine).fillna(0)
    im = ax.imshow(df.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(all_cats_fine)))
    ax.set_yticks(range(len(all_cats_fine)))
    ax.set_xticklabels(all_cats_fine, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(all_cats_fine, fontsize=8)
    n_imm_grp = (imm.obs["phase2_group"] == g).sum()
    ax.set_title(f"{g}   n_immune={n_imm_grp:,}", fontsize=11)
fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, label="nhood enrichment z-score")
fig.suptitle("SCART Phase 2 — Fig 3A(ii): neighborhood enrichment by 14-type immune annotation", fontsize=13, y=1.02)
plt.savefig(FIG / "phase2_fig3A_nhood_immune_per_group.png", dpi=250, bbox_inches="tight")
plt.savefig(FIG / "phase2_fig3A_nhood_immune_per_group.pdf", bbox_inches="tight")
plt.close()
print(f"  [fig] phase2_fig3A_nhood_immune_per_group.png/.pdf")

# ============================================================================
# (C) SCART3_L vs SBRT_L focal delta — L-zone immune cells only
# ============================================================================
print("\n[C] Rim showdown: SCART3_L vs SBRT_L immune enrichment delta")
l_results = {}
for g in ["SBRT", "SCART3"]:
    sub_mask = (imm.obs["phase2_group"] == g) & (imm.obs["phase2_zone"] == "L")
    if sub_mask.sum() < 30:
        print(f"  {g}_L too few cells: {sub_mask.sum()}"); continue
    sub = imm[sub_mask].copy()
    sub.obs[fine_col] = sub.obs[fine_col].astype(str).astype("category")
    sub.obs[fine_col] = sub.obs[fine_col].cat.remove_unused_categories()
    sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=6)
    sq.gr.nhood_enrichment(sub, cluster_key=fine_col, seed=42, n_perms=200, n_jobs=1, show_progress_bar=False)
    z = sub.uns[f"{fine_col}_nhood_enrichment"]["zscore"]
    cats = sub.obs[fine_col].cat.categories.tolist()
    l_results[g] = pd.DataFrame(z, index=cats, columns=cats)
    print(f"  {g}_L: n={sub.n_obs}  clusters={len(cats)}")

if "SBRT" in l_results and "SCART3" in l_results:
    all_l_cats = sorted(set(l_results["SBRT"].index) | set(l_results["SCART3"].index))
    sbrt_z = l_results["SBRT"].reindex(index=all_l_cats, columns=all_l_cats).fillna(0)
    scart3_z = l_results["SCART3"].reindex(index=all_l_cats, columns=all_l_cats).fillna(0)
    delta = scart3_z - sbrt_z

    fig, axes = plt.subplots(1, 3, figsize=(22, 7))
    vmax = max(np.abs(sbrt_z.values).max(), np.abs(scart3_z.values).max(), 3)
    for ax, (title, mat) in zip(axes, [("SBRT_L (rim)", sbrt_z),
                                        ("SCART3_L (rim)", scart3_z),
                                        ("Δ = SCART3_L − SBRT_L", delta)]):
        vm = vmax if "Δ" not in title else np.abs(delta.values).max()
        im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vm, vmax=vm, aspect="auto")
        ax.set_xticks(range(len(all_l_cats)))
        ax.set_yticks(range(len(all_l_cats)))
        ax.set_xticklabels(all_l_cats, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(all_l_cats, fontsize=8)
        ax.set_title(title, fontsize=11)
        plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    fig.suptitle("SCART Phase 2 — Fig 3B(quant): SCART3_L vs SBRT_L neighborhood enrichment", fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(FIG / "phase2_fig3B_scart3L_vs_sbrtL_immune_delta.png", dpi=250, bbox_inches="tight")
    plt.savefig(FIG / "phase2_fig3B_scart3L_vs_sbrtL_immune_delta.pdf", bbox_inches="tight")
    plt.close()
    print(f"  [fig] phase2_fig3B_scart3L_vs_sbrtL_immune_delta.png/.pdf")

    # top enrichment changes
    d = delta.stack().sort_values(ascending=False)
    d = d[~d.index.duplicated()]
    # only keep one side of symmetric pairs
    keep = [(a,b) for a,b in d.index if a < b]
    d = d.loc[keep]
    print("\n  TOP 10 SCART3_L ENRICHED pairs (vs SBRT_L):")
    print(d.head(10).to_string())
    print("\n  TOP 10 SBRT_L ENRICHED pairs (vs SCART3_L):")
    print(d.tail(10).to_string())

# ============================================================================
# save all z-score tables
# ============================================================================
all_z_long = []
for g, df in broad_results.items():
    t = df.stack().reset_index()
    t.columns = ["cluster_a", "cluster_b", "zscore"]
    t["group"] = g; t["annotation"] = "broad"
    all_z_long.append(t)
for g, df in fine_results.items():
    t = df.stack().reset_index()
    t.columns = ["cluster_a", "cluster_b", "zscore"]
    t["group"] = g; t["annotation"] = "fine_immune"
    all_z_long.append(t)
for g, df in l_results.items():
    t = df.stack().reset_index()
    t.columns = ["cluster_a", "cluster_b", "zscore"]
    t["group"] = g + "_L"; t["annotation"] = "fine_immune_L_only"
    all_z_long.append(t)
all_z_df = pd.concat(all_z_long, ignore_index=True)
all_z_df.to_csv(QC / "nhood_enrichment_zscores.csv", index=False)
print(f"\n[done] wrote nhood_enrichment_zscores.csv ({len(all_z_df)} rows)")
