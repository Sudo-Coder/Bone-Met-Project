#!/usr/bin/env python
"""Task 1 downstream — regulon targets, AUCell by group, hypergeometric LA-TAM overlap,
TAM/faction specificity. Reads pySCENIC outputs (regulons.csv, auc_mtx.csv) + cell metadata.
Run AFTER the SLURM pyscenic job:
  envs/mechanism_env/bin/python analysis/03_pyscenic_downstream.py
Writes only to outputs/.
"""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom, mannwhitneyu, spearmanr

ROOT = "/autofs/projects-t3/hussain/scProj"
PROJ = ROOT + "/kidney-cancer"   # project data/outputs (relocated 2026-07-02)
SCE  = f"{PROJ}/outputs/scenic"
TAB  = f"{PROJ}/outputs/tables"; FIG = f"{PROJ}/outputs/figures"
os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)

LATAM = ["TREM2","APOE","APOC1","C1QA","C1QB","C1QC","CD9","FABP5"]
FOCUS_TF = ["NR4A2","NR4A3","NFKB1","RELA","HIF1A","NFE2L2"]

# ---- load regulons ----
from pyscenic.cli.utils import load_signatures
regs = load_signatures(f"{SCE}/regulons.csv")
reg_by_tf = {}
for r in regs:
    tf = r.name.split("(")[0]
    reg_by_tf.setdefault(tf, set()).update(r.genes)
print("n regulons:", len(regs), "| TFs with regulon:", len(reg_by_tf), flush=True)

# target tables for focus TFs
rows = []
for tf in FOCUS_TF:
    genes = sorted(reg_by_tf.get(tf, set()))
    rows.append({"TF": tf, "n_targets": len(genes),
                 "latam_targets": ";".join([g for g in genes if g in LATAM]),
                 "targets": ";".join(genes)})
pd.DataFrame(rows).to_csv(f"{TAB}/scenic_regulon_targets.csv", index=False)
print("wrote scenic_regulon_targets.csv", flush=True)

# ---- hypergeometric overlap of NR4A regulon targets with LA-TAM set ----
universe = pd.read_csv(f"{SCE}/gene_universe.csv")["gene"].tolist()
M = len(universe)                                   # population size
latam_in_univ = [g for g in LATAM if g in universe]
n = len(latam_in_univ)                              # successes in population
hyp_rows = []
for tf in ["NR4A2","NR4A3","NFKB1"]:
    targets = reg_by_tf.get(tf, set()) & set(universe)
    N = len(targets)                                # sample size
    k = len(targets & set(latam_in_univ))           # observed successes
    p = hypergeom.sf(k-1, M, n, N) if N > 0 and n > 0 else np.nan
    hyp_rows.append({"TF": tf, "n_targets_in_universe": N, "LA_TAM_universe": n,
                     "overlap_k": k, "overlap_genes": ";".join(sorted(targets & set(latam_in_univ))),
                     "expected": (n*N/M) if M else np.nan, "hypergeom_p": p})
pd.DataFrame(hyp_rows).to_csv(f"{TAB}/scenic_nr4a_latam_hypergeom.csv", index=False)
print("wrote scenic_nr4a_latam_hypergeom.csv:\n", pd.DataFrame(hyp_rows), flush=True)

# ---- AUCell by group ----
auc = pd.read_csv(f"{SCE}/auc_mtx.csv", index_col=0)
auc.columns = [c.split("(")[0] for c in auc.columns]   # strip (+)
meta = pd.read_parquet(f"{SCE}/cell_metadata.parquet")
common = auc.index.intersection(meta.index)
auc = auc.loc[common]; meta = meta.loc[common]
print("AUCell matrix:", auc.shape, "| cells matched:", len(common), flush=True)

present_tf = [tf for tf in FOCUS_TF if tf in auc.columns]
df = meta.join(auc[present_tf])

# mean AUCell by final_label and by condition
df.groupby("final_label")[present_tf].mean().to_csv(f"{TAB}/scenic_aucell_by_celltype.csv")
df.groupby("condition")[present_tf].mean().to_csv(f"{TAB}/scenic_aucell_by_condition.csv")
df.groupby(["final_label","condition"])[present_tf].mean().to_csv(f"{TAB}/scenic_aucell_by_group.csv")

# ---- specificity: TAM(Tumor) vs rest, per focus TF ----
is_tam_tumor = (df["final_label"]=="TAM") & (df["condition"]=="Tumor")
spec = []
for tf in present_tf:
    a = df.loc[is_tam_tumor, tf].values
    b = df.loc[~is_tam_tumor, tf].values
    if len(a) > 5 and len(b) > 5:
        U, p = mannwhitneyu(a, b, alternative="greater")
        spec.append({"TF": tf, "mean_TAM_tumor": float(np.mean(a)), "mean_rest": float(np.mean(b)),
                     "log2FC": float(np.log2((np.mean(a)+1e-9)/(np.mean(b)+1e-9))),
                     "mannwhitney_U": float(U), "p_greater": float(p), "n_tam_tumor": int(len(a))})
pd.DataFrame(spec).to_csv(f"{TAB}/scenic_tam_tumor_specificity.csv", index=False)
print("wrote scenic_tam_tumor_specificity.csv:\n", pd.DataFrame(spec), flush=True)

# ---- correlation of NR4A AUCell vs CD9/TREM2 expression (per cell) ----
cor = []
for tf in [t for t in ["NR4A2","NR4A3"] if t in present_tf]:
    for g in ["CD9","TREM2","GPNMB","APOE"]:
        col = f"expr_{g}"
        if col in df.columns:
            rho, p = spearmanr(df[tf], df[col])
            cor.append({"TF_AUCell": tf, "gene_expr": g, "spearman_rho": rho, "p": p})
pd.DataFrame(cor).to_csv(f"{TAB}/scenic_nr4a_vs_latam_expr_corr.csv", index=False)

# ---- figures ----
# heatmap: mean AUCell (focus TFs) by celltype x condition
piv = df.groupby(["final_label","condition"])[present_tf].mean()
fig, axes = plt.subplots(1, len(present_tf), figsize=(3.2*len(present_tf), 4), squeeze=False)
for i, tf in enumerate(present_tf):
    m = piv[tf].unstack("condition")
    sns.heatmap(m, ax=axes[0][i], cmap="magma", cbar_kws={"label":"mean AUCell"})
    axes[0][i].set_title(tf); axes[0][i].set_xlabel(""); axes[0][i].set_ylabel("")
plt.tight_layout(); plt.savefig(f"{FIG}/scenic_aucell_heatmap.png", dpi=150); plt.close()

# violins of NR4A/NFKB AUCell by celltype (Tumor faction only where relevant)
for tf in present_tf:
    plt.figure(figsize=(7,3.5))
    order = ["Mono1","Mono2","Mono3","Monocyte Pro","TIM","TAM","mDC","Osteoclasts"]
    order = [o for o in order if o in df["final_label"].unique()]
    sns.violinplot(data=df, x="final_label", y=tf, order=order, cut=0, scale="width")
    plt.xticks(rotation=45, ha="right"); plt.title(f"{tf} regulon AUCell by cell type")
    plt.tight_layout(); plt.savefig(f"{FIG}/scenic_aucell_violin_{tf}.png", dpi=150); plt.close()

print("TASK1 DOWNSTREAM DONE", flush=True)
