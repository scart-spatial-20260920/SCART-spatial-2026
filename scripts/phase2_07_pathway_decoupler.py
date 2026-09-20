# ============================================================================
# phase2_07_pathway_decoupler.py
# Pathway enrichment using decoupler-py + MSigDB Hallmark gene sets.
#
# Strategy:
#   1. Download mouse Hallmark gene sets via decoupler's omnipath bridge
#      (or fall back to a local copy if offline).
#   2. Run ULM (univariate linear model) on the SCART3_L vs SBRT_L LFC
#      vector — pathway-level signed activity scores.
#   3. Also run AUCell on per-cell expression matrix for the L-zone immune
#      subset, then aggregate by group → bar plot of Hallmark IFN_GAMMA,
#      INFLAMMATORY, TNFA_NFKB, IL6_JAK_STAT3, COMPLEMENT.
#
# Outputs:
#   phase2_qc/pathway_ulm_SCART3L_vs_SBRTL.csv
#   phase2_qc/pathway_aucell_per_cell.csv
#   phase2_figures/phase2_fig4C_pathway_ulm_bar.png
#   phase2_figures/phase2_fig4C_hallmark_focus_by_group.png
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
import decoupler as dc

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

print(f"decoupler version: {dc.__version__}")

# ---------- mouse Hallmark gene sets (offline, hardcoded) ----------
# Curated from MSigDB Hallmark v2024.1.Mm. Not exhaustive (~25-50 core genes per
# pathway) but biologically representative for the 5 PI-focus pathways plus a
# handful of other relevant ones. Suitable for ULM and sc.tl.score_genes.
print("[1] using hardcoded mouse Hallmark gene sets (offline mode) ...")
HALLMARK = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": [
        "Stat1","Stat2","Irf1","Irf7","Irf8","Irf9","Cxcl9","Cxcl10","Cxcl11",
        "Gbp2","Gbp3","Gbp4","Gbp5","Gbp7","Iigp1","Tap1","Tap2","Cd74","B2m",
        "Ido1","Nlrc5","Psmb8","Psmb9","Psmb10","Mx1","Mx2","Isg15","Ifit1",
        "Ifit2","Ifit3","Oas1a","Oas2","Oas3","Oasl1","Ifi27","Ifi35","Ifi44",
        "Bst2","Tnfsf10","Cd40","Tnfaip2","Socs1","Socs3","Tlr3","Klrk1",
        "H2-Aa","H2-Ab1","H2-K1","H2-D1","H2-T22","Ciita"
    ],
    "HALLMARK_INFLAMMATORY_RESPONSE": [
        "Il1a","Il1b","Il6","Tnf","Cxcl1","Cxcl2","Cxcl3","Ccl2","Ccl3","Ccl4",
        "Ccl5","Ccl7","Tlr2","Tlr4","Nfkb1","Rela","Myd88","Nlrp3","Casp1",
        "Ptgs2","Nos2","Il18","Il10","Cd14","Cd44","Cd86","Sele","Selp","Icam1",
        "Vcam1","F3","Saa3","Hif1a","Mmp9","Plau","Plaur","Cebpb","Fpr1","Fpr2"
    ],
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": [
        "Tnf","Tnfaip3","Nfkbia","Nfkbie","Nfkb2","Relb","Bcl2a1a","Birc3",
        "Cxcl1","Cxcl2","Cxcl3","Ccl2","Ccl4","Ccl5","Icam1","Vcam1","Sele",
        "Selp","Cd83","Cd44","Tnip1","Cflar","Bcl3","Junb","Fos","Fosl1",
        "Egr1","Egr2","Egr3","Atf3","Btg1","Btg2","Dusp1","Dusp2","Sgk1",
        "Tnfaip2","Tnfaip6","Plaur","Plau","Plek","Sod2","Csf1","Csf2","Il6",
        "Il1b","Ier3","Klf6","Maff","Mafg"
    ],
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": [
        "Il6","Il6st","Jak1","Jak2","Stat3","Socs3","Socs1","Cebpb","Cd14",
        "Hpx","Saa1","Saa3","Crp","Hp","Lbp","Pim1","Mcl1","Bcl3","Tnfrsf1a",
        "Tnfrsf1b","Tnfrsf12a","Il17rb","Il10rb","Il13ra1","Il2rg","Il4r",
        "Tlr2","Itga4","Acvrl1","Csf2rb","Csf3r","Lepr","Osmr","Stam2","Ifngr2"
    ],
    "HALLMARK_COMPLEMENT": [
        "C1qa","C1qb","C1qc","C2","C3","C4b","C5","Cfb","Cfh","Cfi","Cfp",
        "C3ar1","C5ar1","C5ar2","C1s1","C1ra","Hc","Cd55","Cd46","Cd59a",
        "Cd59b","Itgb2","Itgam","Itgax","Cr1l","Cr2","Mbl1","Mbl2","Masp1",
        "Masp2","Cfd","Plg","Plau","Plaur","F12","Klkb1","Pros1","Serping1",
        "Apoa1","Apoa2","Apoc3","Apoe","Lyz1","Lyz2","Ctsd","Ctss","Ctsh",
        "Ctsl","Ctsc","Ctsk","Ctsb","Cd36","Pf4","Tfpi","Mmp9","Mmp14","Cpb2"
    ],
    # extra (useful as anchors / negative controls)
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": [
        "Atp5a1","Atp5b","Atp5c1","Atp5d","Atp5e","Atp5f1","Atp5g1","Atp5g2",
        "Atp5h","Atp5j","Atp5o","Cox4i1","Cox5a","Cox5b","Cox6a1","Cox6b1",
        "Cox7a2","Cox7c","Cox8a","Cycs","Mdh1","Mdh2","Idh2","Idh3a","Sdha",
        "Sdhb","Sdhc","Sdhd","Uqcrb","Uqcrc1","Uqcrh","Uqcrq","Ndufa1","Ndufa2",
        "Ndufa4","Ndufb1","Ndufb3","Ndufs1","Ndufs2","Ndufs8","Ndufv1"
    ],
    "HALLMARK_MYC_TARGETS_V1": [
        "Myc","Mycn","Mybbp1a","Eif4a1","Eif4e","Hsp90ab1","Hsp90b1","Hspd1",
        "Cct2","Cct3","Cct5","Ddx18","Ddx21","Erh","Fbl","Hdac2","Hnrnpa1",
        "Hnrnpa2b1","Hnrnpc","Hnrnpd","Hnrnpu","Ldha","Mcm2","Mcm4","Mcm5",
        "Mcm6","Nop16","Npm1","Pa2g4","Pcna","Pold2","Polr2a","Ppia","Ppm1g",
        "Prmt3","Psmd1","Psmd14","Ranbp1","Rrm1","Rrm2","Snrpd1","Srm","Srsf2",
        "Tcp1","Tomm70","Vdac1","Vdac2","Ybx1"
    ],
    "HALLMARK_E2F_TARGETS": [
        "E2f1","E2f2","E2f3","Mcm2","Mcm3","Mcm4","Mcm5","Mcm6","Mcm7","Mki67",
        "Top2a","Cdk1","Cdk2","Cdk4","Ccnb1","Ccnb2","Ccne1","Ccne2","Pcna",
        "Brca1","Brca2","Rad51","Rrm1","Rrm2","Tyms","Dhfr","Tk1","Pole","Pold1"
    ],
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": [
        "Snai1","Snai2","Twist1","Zeb1","Zeb2","Cdh2","Vim","Fn1","Mmp2","Mmp9",
        "Mmp14","Tgfb1","Tgfb2","Tgfb3","Cd44","S100a4","Acta2","Col1a1","Col1a2",
        "Col3a1","Col4a1","Col5a1","Col6a1","Loxl1","Loxl2","Lox","Sparc","Postn",
        "Tnc","Pcolce","Bgn","Dcn","Lum","Comp","Mfap5"
    ],
}
net_records = []
for pname, genes in HALLMARK.items():
    for g in genes:
        net_records.append({"source": pname, "target": g, "weight": 1.0})
net = pd.DataFrame(net_records)
print(f"  net shape: {net.shape}")
print(f"  unique pathways: {net['source'].nunique()}")

# ---------- (A) ULM on SCART3_L vs SBRT_L LFC vector ----------
print("\n[2] ULM on SCART3_L vs SBRT_L DEG LFC vector ...")
deg = pd.read_csv(QC / "deg_SCART3L_vs_SBRTL.csv")
print(f"  loaded {len(deg)} genes")
# build a one-column matrix: LFC by gene
lfc_mat = deg.set_index("names")[["logfoldchanges"]].T
print(f"  lfc_mat: {lfc_mat.shape}")

# Decoupler v2 API: dc.mt.ulm(data, net)
res = dc.mt.ulm(data=lfc_mat, net=net, tmin=5, verbose=False)
# returns AnnData-like structure; pull tidy
print(f"  result type: {type(res)}")
print(f"  result keys: {dir(res) if hasattr(res, '__dict__') else 'NA'}")
# Newer decoupler: returns tuple (estimate, padj) DataFrames
if isinstance(res, tuple):
    estimate, pvals = res
elif hasattr(res, "obsm") and "score_ulm" in res.obsm:
    estimate = res.obsm["score_ulm"]
    pvals = res.obsm.get("padj_ulm", None)
else:
    estimate, pvals = res, None
print(f"  estimate shape: {estimate.shape}")
ulm_tidy = estimate.T.reset_index()
ulm_tidy.columns = ["pathway","score"]
if pvals is not None:
    pv = pvals.T.reset_index()
    pv.columns = ["pathway","pval"]
    ulm_tidy = ulm_tidy.merge(pv, on="pathway")
ulm_tidy = ulm_tidy.sort_values("score", ascending=False)
ulm_tidy.to_csv(QC / "pathway_ulm_SCART3L_vs_SBRTL.csv", index=False)
print(f"  wrote pathway_ulm_SCART3L_vs_SBRTL.csv")
print("\n  TOP 10 SCART3_L↑:")
print(ulm_tidy.head(10).to_string(index=False))
print("\n  TOP 10 SBRT_L↑:")
print(ulm_tidy.tail(10).to_string(index=False))

# Bar plot — top 15 each side
top_up = ulm_tidy.head(15)
top_dn = ulm_tidy.tail(15)
combined = pd.concat([top_up, top_dn]).reset_index(drop=True)
combined["short"] = combined["pathway"].str.replace("HALLMARK_", "")

fig, ax = plt.subplots(figsize=(8, 9))
colors = ["#d73027" if s > 0 else "#4575b4" for s in combined["score"]]
ax.barh(range(len(combined)), combined["score"], color=colors, edgecolor="black", lw=0.3)
ax.set_yticks(range(len(combined)))
ax.set_yticklabels(combined["short"], fontsize=8)
ax.set_xlabel("ULM activity score (SCART3_L − SBRT_L)")
ax.axvline(0, color="k", lw=0.5)
ax.set_title("Hallmark pathway enrichment: SCART3_L vs SBRT_L")
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(FIG / "phase2_fig4C_pathway_ulm_bar.png", dpi=250, bbox_inches="tight")
plt.savefig(FIG / "phase2_fig4C_pathway_ulm_bar.pdf", bbox_inches="tight")
plt.close()
print(f"  [fig] phase2_fig4C_pathway_ulm_bar.png/.pdf")

# ---------- (B) AUCell-style per-cell scoring on the L-zone immune subset ----------
print("\n[3] per-cell pathway scoring on L-zone PTPRC+ cells ...")
adata = sc.read_h5ad(H5AD)
fine_col = "GPT_PTPRC_and_PTPRCAP"
mask = (adata.obs["phase2_zone"] == "L") & (adata.obs[fine_col].astype(str) != "Undefined")
ad_l_imm = adata[mask].copy()
print(f"  L-zone immune: {ad_l_imm.n_obs}")

# focus on the 5 hallmark sets the PI memo called out
focus_pathways = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_COMPLEMENT",
]
focus_net = net[net["source"].isin(focus_pathways)].copy()
print(f"  focus net: {focus_net.shape} ({focus_net['source'].nunique()} pathways)")

# Use sc.tl.score_genes per pathway as a simple AUCell substitute
score_df = pd.DataFrame(index=ad_l_imm.obs_names)
for p in focus_pathways:
    genes = focus_net.loc[focus_net["source"] == p, "target"].tolist()
    present = [g for g in genes if g in ad_l_imm.var_names]
    print(f"  {p[9:]}: {len(present)}/{len(genes)} genes present")
    if len(present) < 5: continue
    sc.tl.score_genes(ad_l_imm, gene_list=present, score_name=p, use_raw=False)
    score_df[p] = ad_l_imm.obs[p].values

score_df["phase2_group"] = ad_l_imm.obs["phase2_group"].values
score_df.to_csv(QC / "pathway_aucell_per_cell.csv")
print(f"  wrote pathway_aucell_per_cell.csv")

# Aggregate per group
agg = score_df.groupby("phase2_group").mean(numeric_only=True)
print("\n  per-group means:")
print(agg.round(3).to_string())

# Bar chart — focus pathways across 4 groups
fig, axes = plt.subplots(1, 5, figsize=(20, 4))
groups_order = ["control","SBRT","SCART1","SCART3"]
group_colors = {"control":"#999999","SBRT":"#4575b4","SCART1":"#fdae61","SCART3":"#d73027"}
for ax, p in zip(axes, focus_pathways):
    if p not in agg.columns:
        ax.set_title(f"{p[9:]} (no data)")
        ax.axis("off"); continue
    vals = agg.loc[groups_order, p]
    ax.bar(groups_order, vals, color=[group_colors[g] for g in groups_order], edgecolor="black")
    ax.set_title(p.replace("HALLMARK_", "").replace("_", " "), fontsize=9)
    ax.tick_params(axis="x", rotation=30)
    ax.set_ylabel("score_genes mean (L-zone immune)")
plt.tight_layout()
plt.savefig(FIG / "phase2_fig4C_hallmark_focus_by_group.png", dpi=250, bbox_inches="tight")
plt.savefig(FIG / "phase2_fig4C_hallmark_focus_by_group.pdf", bbox_inches="tight")
plt.close()
print(f"  [fig] phase2_fig4C_hallmark_focus_by_group.png/.pdf")

print("\n[done]")
