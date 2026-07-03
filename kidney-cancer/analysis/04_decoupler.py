#!/usr/bin/env python
"""Task 3 — decoupleR PROGENy (pathways) + CollecTRI (TF activity) per cell.
API: decoupler 2.1.4 (dc.op.progeny/collectri nets; dc.mt.mlm/ulm methods).
Run: envs/mechanism_env/bin/python analysis/04_decoupler.py
Reads myeloid h5ad (read-only) + optional pseudotime parquet; writes only outputs/.
"""
import os, numpy as np, pandas as pd, scanpy as sc, decoupler as dc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns

ROOT="/autofs/projects-t3/hussain/scProj"
PROJ = ROOT + "/kidney-cancer"   # project data/outputs (relocated 2026-07-02)
H5AD=f"{PROJ}/Cleaned_Data/myeloid_FINAL_labels.h5ad"
TAB=f"{PROJ}/outputs/tables"; FIG=f"{PROJ}/outputs/figures"; DEC=f"{PROJ}/outputs/decoupler"
for d in (TAB,FIG,DEC): os.makedirs(d, exist_ok=True)
PT=f"{PROJ}/outputs/cellrank/pseudotime.parquet"   # produced by task 4 (optional)

FOCUS_PATH=["NFkB","Hypoxia","JAK-STAT","TNFa","TGFb","MAPK","p53","VEGF"]
FOCUS_TF=["NR4A2","NR4A3","NFKB1","RELA","HIF1A","NFE2L2"]

print("loading", flush=True)
adata=sc.read_h5ad(H5AD)
adata.X=adata.layers["log_norm"].copy()          # log-normalized for activity inference
adata.obs["final_label"]=adata.obs["final_label"].astype(str)
adata.obs["condition"]=adata.obs["condition"].astype(str)

def run(net, method, tag):
    fn=dict(mlm=dc.mt.mlm, ulm=dc.mt.ulm)[method]
    fn(data=adata, net=net, tmin=5, verbose=True)
    skey=[k for k in adata.obsm if k.startswith("score_")][-1]
    scores=adata.obsm[skey].copy()
    del adata.obsm[skey]
    pkey=[k for k in adata.obsm if k.startswith("padj_")]
    if pkey: del adata.obsm[pkey[-1]]
    scores.index=adata.obs_names
    scores.to_parquet(f"{DEC}/{tag}_percell.parquet")
    scores.to_csv(f"{TAB}/decoupler_{tag}_percell.csv")
    print(f"{tag}: {scores.shape}, sources:", list(scores.columns)[:12], flush=True)
    return scores

progeny=dc.op.progeny(organism="human", top=500)
collectri=dc.op.collectri(organism="human")
print("progeny net:", progeny.shape, "| collectri net:", collectri.shape, flush=True)

path_act=run(progeny, "mlm", "progeny")
tf_act  =run(collectri, "ulm", "collectri")

# ---- grouped summaries (final_label x condition) ----
meta=adata.obs[["final_label","condition"]].copy()
def grouped(scores, cols, tag):
    cols=[c for c in cols if c in scores.columns]
    df=meta.join(scores[cols])
    df.groupby("final_label")[cols].mean().to_csv(f"{TAB}/decoupler_{tag}_by_celltype.csv")
    df.groupby("condition")[cols].mean().to_csv(f"{TAB}/decoupler_{tag}_by_condition.csv")
    df.groupby(["final_label","condition"])[cols].mean().to_csv(f"{TAB}/decoupler_{tag}_by_group.csv")
    return df, cols
pdf,pcols=grouped(path_act, FOCUS_PATH, "progeny")
tdf,tcols=grouped(tf_act,   FOCUS_TF,   "collectri")

# ---- figures: heatmaps celltype x condition ----
def heatmaps(df, cols, tag):
    piv=df.groupby(["final_label","condition"])[cols].mean()
    fig,axes=plt.subplots(1,len(cols),figsize=(2.9*len(cols),4),squeeze=False)
    for i,c in enumerate(cols):
        m=piv[c].unstack("condition")
        sns.heatmap(m,ax=axes[0][i],cmap="coolwarm",center=0,cbar_kws={"label":"activity"})
        axes[0][i].set_title(c); axes[0][i].set_xlabel(""); axes[0][i].set_ylabel("")
    plt.tight_layout(); plt.savefig(f"{FIG}/decoupler_{tag}_heatmap.png",dpi=150); plt.close()
heatmaps(pdf,pcols,"progeny"); heatmaps(tdf,tcols,"collectri")

# ---- boxplots by faction for key signals (TAM cells only, across factions) ----
for scores,cols,tag in [(path_act,["NFkB","Hypoxia","JAK-STAT","TNFa"],"progeny"),
                        (tf_act,["NR4A2","NFKB1","HIF1A","RELA"],"collectri")]:
    cols=[c for c in cols if c in scores.columns]
    df=meta.join(scores[cols])
    fig,axes=plt.subplots(1,len(cols),figsize=(3.2*len(cols),3.6),squeeze=False)
    order=["Benign","Distal","Involved","Tumor"]
    for i,c in enumerate(cols):
        sns.boxplot(data=df,x="condition",y=c,order=order,ax=axes[0][i],showfliers=False)
        axes[0][i].set_title(c); axes[0][i].set_xlabel(""); axes[0][i].tick_params(axis='x',rotation=45)
    plt.tight_layout(); plt.savefig(f"{FIG}/decoupler_{tag}_byfaction.png",dpi=150); plt.close()

# ---- along pseudotime (if task 4 produced it) ----
if os.path.exists(PT):
    pt=pd.read_parquet(PT)
    col=[c for c in pt.columns if "pseudotime" in c.lower()]
    if col:
        pt=pt[col[0]].rename("pseudotime")
        for scores,cols,tag in [(path_act,["NFkB","Hypoxia","JAK-STAT","TNFa"],"progeny"),
                                (tf_act,["NR4A2","NR4A3","NFKB1","HIF1A"],"collectri")]:
            cols=[c for c in cols if c in scores.columns]
            j=scores[cols].join(pt,how="inner").dropna(subset=["pseudotime"]).sort_values("pseudotime")
            plt.figure(figsize=(7,4))
            for c in cols:
                y=j[c].rolling(200,min_periods=20,center=True).mean()
                plt.plot(j["pseudotime"],y,label=c)
            plt.legend(); plt.xlabel("pseudotime"); plt.ylabel("activity (rolling mean)")
            plt.title(f"{tag} activity along Mono1->TAM pseudotime")
            plt.tight_layout(); plt.savefig(f"{FIG}/decoupler_{tag}_along_pseudotime.png",dpi=150); plt.close()
    print("pseudotime overlays written", flush=True)
else:
    print("NOTE: pseudotime parquet not found yet; run task 4 then re-run for pseudotime overlays.", flush=True)

print("TASK3 DECOUPLER DONE", flush=True)
