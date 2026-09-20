# ============================================================================
# phase2_06_deg_rim.py
# Differential gene expression: SCART3_L vs SBRT_L (and a few other contrasts)
#
# Why this comparison: same batch (_D1_D31), both rim, directly addresses
# "is SCART3 rim more antigen-presenting than SBRT rim".
#
# Outputs:
#   phase2_qc/deg_SCART3L_vs_SBRTL.csv
#   phase2_qc/deg_SCART3L_vs_controlL.csv
#   phase2_qc/deg_SBRTL_vs_controlL.csv
#   phase2_figures/phase2_fig4B_volcano_SCART3L_vs_SBRTL.png
#   phase2_figures/phase2_fig4B_top_genes_heatmap.png
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

# unbuffered prints
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
print(f"  full: {adata.n_obs:,} × {adata.n_vars:,}")

# Subset to L-zone IMMUNE (PTPRC+) cells across all groups.
# Reason: Phase 1 was computed on PTPRC+ subset; running DEG on all cells dilutes
# immune signal with stroma and produces stromal contamination artifacts.
fine_col = "GPT_PTPRC_and_PTPRCAP"
imm_mask = adata.obs[fine_col].astype(str) != "Undefined"
l_mask = (adata.obs["phase2_zone"] == "L") & imm_mask
print(f"  L-zone IMMUNE cells (PTPRC+): {l_mask.sum():,}")
ad_l = adata[l_mask].copy()
ad_l.obs["phase2_group"] = ad_l.obs["phase2_group"].astype(str).astype("category")
ad_l.obs["phase2_group"] = ad_l.obs["phase2_group"].cat.remove_unused_categories()
print(f"  L-zone group counts:")
print(ad_l.obs["phase2_group"].value_counts())

# scanpy's rank_genes_groups expects log-normalized X. We loaded h5ad with X = logcounts (Seurat 'data' layer)
# Verify
print(f"  X dtype: {ad_l.X.dtype}, max: {ad_l.X.max():.2f}, min: {ad_l.X.min():.2f}")

# Pre-filter to genes expressed in ≥3% of cells in at least one of the groups
# Reason: with 26,757 sparse genes and small n, FDR correction wipes everything.
# Filtering to expressed-enough genes recovers power for the immune signal.
print("[filter] keeping genes with ≥3% expression in any group...")
keep = np.zeros(ad_l.n_vars, dtype=bool)
for g in ad_l.obs["phase2_group"].unique():
    sub = ad_l[ad_l.obs["phase2_group"] == g]
    if sub.n_obs == 0: continue
    Xs = sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X
    pct = (Xs > 0).mean(axis=0)
    keep |= (pct >= 0.03)
print(f"  kept {keep.sum():,} / {ad_l.n_vars:,} genes")
ad_l = ad_l[:, keep].copy()

# ============================================================================
# DEG: SCART3 vs SBRT in L-zone
# ============================================================================
print("\n[DEG] SCART3_L vs SBRT_L (Wilcoxon)")
ad_pair = ad_l[ad_l.obs["phase2_group"].isin(["SCART3", "SBRT"])].copy()
ad_pair.obs["phase2_group"] = ad_pair.obs["phase2_group"].astype(str).astype("category")
ad_pair.obs["phase2_group"] = ad_pair.obs["phase2_group"].cat.remove_unused_categories()
print(f"  comparing: {dict(ad_pair.obs['phase2_group'].value_counts())}")

sc.tl.rank_genes_groups(ad_pair, "phase2_group", reference="SBRT",
                         method="wilcoxon", pts=True, use_raw=False)
res = sc.get.rank_genes_groups_df(ad_pair, group="SCART3")
res = res.sort_values("logfoldchanges", ascending=False)
res.to_csv(QC / "deg_SCART3L_vs_SBRTL.csv", index=False)
print(f"  wrote deg_SCART3L_vs_SBRTL.csv ({len(res)} genes)")
print("  TOP 15 SCART3↑:")
print(res.head(15)[[c for c in ["names","logfoldchanges","pvals_adj","pct_nz_group","pct_nz_reference"] if c in res.columns]].to_string())
print("\n  TOP 15 SBRT↑ (SCART3↓):")
print(res.tail(15)[[c for c in ["names","logfoldchanges","pvals_adj","pct_nz_group","pct_nz_reference"] if c in res.columns]].to_string())

# ============================================================================
# DEG: SCART3 vs control in L-zone
# ============================================================================
print("\n[DEG] SCART3_L vs control_L (Wilcoxon)  — note: cross-batch")
ad_pair2 = ad_l[ad_l.obs["phase2_group"].isin(["SCART3","control"])].copy()
ad_pair2.obs["phase2_group"] = ad_pair2.obs["phase2_group"].astype(str).astype("category")
sc.tl.rank_genes_groups(ad_pair2, "phase2_group", reference="control",
                         method="wilcoxon", pts=True, use_raw=False)
res2 = sc.get.rank_genes_groups_df(ad_pair2, group="SCART3")
res2 = res2.sort_values("logfoldchanges", ascending=False)
res2.to_csv(QC / "deg_SCART3L_vs_controlL.csv", index=False)
print(f"  wrote deg_SCART3L_vs_controlL.csv ({len(res2)} genes)")

# ============================================================================
# DEG: SBRT vs control in L-zone
# ============================================================================
print("\n[DEG] SBRT_L vs control_L (Wilcoxon)  — note: cross-batch")
ad_pair3 = ad_l[ad_l.obs["phase2_group"].isin(["SBRT","control"])].copy()
ad_pair3.obs["phase2_group"] = ad_pair3.obs["phase2_group"].astype(str).astype("category")
sc.tl.rank_genes_groups(ad_pair3, "phase2_group", reference="control",
                         method="wilcoxon", pts=True, use_raw=False)
res3 = sc.get.rank_genes_groups_df(ad_pair3, group="SBRT")
res3 = res3.sort_values("logfoldchanges", ascending=False)
res3.to_csv(QC / "deg_SBRTL_vs_controlL.csv", index=False)
print(f"  wrote deg_SBRTL_vs_controlL.csv ({len(res3)} genes)")

# ============================================================================
# Volcano plot — SCART3_L vs SBRT_L
# ============================================================================
print("\n[volcano] plotting ...")
v = res.copy()
v["minus_log10_padj"] = -np.log10(v["pvals_adj"].clip(lower=1e-300))
# truncate extreme values
v.loc[v["minus_log10_padj"] > 50, "minus_log10_padj"] = 50

sig = (v["pvals_adj"] < 0.05) & (v["logfoldchanges"].abs() > 0.5)
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(v.loc[~sig, "logfoldchanges"], v.loc[~sig, "minus_log10_padj"],
           s=4, c="#cccccc", alpha=0.5, rasterized=True)
ax.scatter(v.loc[sig & (v["logfoldchanges"]>0), "logfoldchanges"],
           v.loc[sig & (v["logfoldchanges"]>0), "minus_log10_padj"],
           s=8, c="#d73027", alpha=0.85, label="↑ in SCART3_L", rasterized=True)
ax.scatter(v.loc[sig & (v["logfoldchanges"]<0), "logfoldchanges"],
           v.loc[sig & (v["logfoldchanges"]<0), "minus_log10_padj"],
           s=8, c="#4575b4", alpha=0.85, label="↑ in SBRT_L", rasterized=True)

# label top 10 each side
top_up = v.loc[sig & (v["logfoldchanges"]>0)].nlargest(12, "logfoldchanges")
top_dn = v.loc[sig & (v["logfoldchanges"]<0)].nsmallest(12, "logfoldchanges")
for _, row in pd.concat([top_up, top_dn]).iterrows():
    ax.annotate(row["names"], (row["logfoldchanges"], row["minus_log10_padj"]),
                fontsize=7, alpha=0.9)

ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.5)
ax.axvline(0, color="grey", ls="--", lw=0.5)
ax.set_xlabel("log2 fold change (SCART3_L / SBRT_L)")
ax.set_ylabel("-log10(adj p-value)")
ax.set_title(f"DEG: SCART3_L vs SBRT_L  (n_SCART3={(ad_pair.obs['phase2_group']=='SCART3').sum()}, n_SBRT={(ad_pair.obs['phase2_group']=='SBRT').sum()})")
ax.legend(loc="upper left", fontsize=9)
plt.tight_layout()
plt.savefig(FIG / "phase2_fig4B_volcano_SCART3L_vs_SBRTL.png", dpi=250, bbox_inches="tight")
plt.savefig(FIG / "phase2_fig4B_volcano_SCART3L_vs_SBRTL.pdf", bbox_inches="tight")
plt.close()
print(f"  [fig] phase2_fig4B_volcano_SCART3L_vs_SBRTL.png/.pdf")

# ============================================================================
# Top-gene heatmap — top 30 each direction across all 4 L-zone groups
# ============================================================================
print("\n[heatmap] top genes across all L-zone groups ...")
top_genes = pd.concat([res.head(20)["names"], res.tail(20)["names"]]).tolist()
present = [g for g in top_genes if g in ad_l.var_names]
print(f"  using {len(present)} genes")

# mean log-expression per group
mean_df = pd.DataFrame(index=present)
for g in ["control","SBRT","SCART1","SCART3"]:
    sub = ad_l[ad_l.obs["phase2_group"] == g, present]
    if sub.n_obs == 0:
        mean_df[g] = np.nan; continue
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X
    mean_df[g] = X.mean(axis=0)
# z-score per gene across groups
mean_z = mean_df.sub(mean_df.mean(axis=1), axis=0).div(mean_df.std(axis=1).replace(0, 1), axis=0)

fig, ax = plt.subplots(figsize=(5, max(8, len(present)*0.18)))
im = ax.imshow(mean_z.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
ax.set_xticks(range(mean_z.shape[1]))
ax.set_xticklabels(mean_z.columns, rotation=45, ha="right")
ax.set_yticks(range(mean_z.shape[0]))
ax.set_yticklabels(mean_z.index, fontsize=7)
ax.set_title("Top DEG (SCART3_L vs SBRT_L)\nmean log-expr z-score across L-zone groups")
plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="z-score")
plt.tight_layout()
plt.savefig(FIG / "phase2_fig4B_top_genes_heatmap.png", dpi=250, bbox_inches="tight")
plt.savefig(FIG / "phase2_fig4B_top_genes_heatmap.pdf", bbox_inches="tight")
plt.close()
print(f"  [fig] phase2_fig4B_top_genes_heatmap.png/.pdf")

# ============================================================================
# Phase 1 immune panel check
# ============================================================================
panel = ["H2-Ab1","H2-Aa","Ccr7","Cd86","Cd74",
         "Isg15","Ifit1","Ifit3","Cxcl10","Cxcl9","Ifng",
         "Foxp3","Il10","Tgfb1","Pdcd1","Lag3",
         "Nos2","Il1b","Arg1","Mrc1"]
panel_present = [g for g in panel if g in ad_l.var_names]
print(f"\n[panel] {len(panel_present)}/{len(panel)} of Phase 1 immune panel found in data")
panel_missing = set(panel) - set(panel_present)
if panel_missing:
    print(f"  missing: {panel_missing}")

panel_means = pd.DataFrame(index=panel_present)
for g in ["control","SBRT","SCART1","SCART3"]:
    sub = ad_l[ad_l.obs["phase2_group"] == g, panel_present]
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X
    panel_means[f"{g}_L"] = X.mean(axis=0)
panel_means.to_csv(QC / "phase1_panel_means_in_L_zone.csv")
print(f"  wrote phase1_panel_means_in_L_zone.csv")
print(panel_means.round(3).to_string())

print("\n[done]")
