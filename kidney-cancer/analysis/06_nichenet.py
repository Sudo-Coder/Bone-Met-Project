#!/usr/bin/env python
"""Task 2 — NicheNet (Python reimplementation) niche -> TAM.
Senders = non-myeloid niche cells (Tumor, MSC-1/2/3, Peri-1/2/3, Endothelial, osteoblasts,
osteoclasts) from integrated.h5ad, labelled via Cleaned_Data/cell-annotations.csv.
Receiver = TAM (myeloid_FINAL_labels.h5ad, final_label==TAM).
Geneset of interest = LA-TAM program (TAM vs Mono1 up-regulated).
Ligand activity = Pearson corr of each ligand's target regulatory-potential vector with the
binary geneset over the receiver background (NicheNet's primary metric) + AUPR.
Run on a compute node (loads integrated.h5ad in memory):
  envs/mechanism_env/bin/python analysis/06_nichenet.py
Writes only outputs/.
"""
import os, numpy as np, pandas as pd, scanpy as sc
import scipy.sparse as sp
from scipy.stats import pearsonr
from sklearn.metrics import average_precision_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns

ROOT="/autofs/projects-t3/hussain/scProj"
PROJ = ROOT + "/kidney-cancer"   # project data/outputs (relocated 2026-07-02)
NN=f"{ROOT}/resources/nichenet"; TAB=f"{PROJ}/outputs/tables"; FIG=f"{PROJ}/outputs/figures"
OUT=f"{PROJ}/outputs/nichenet"
for d in (TAB,FIG,OUT): os.makedirs(d, exist_ok=True)

LATAM=["TREM2","APOE","APOC1","C1QA","C1QB","C1QC","CD9","FABP5"]
SENDER_TYPES=["Tumor","MSC-1","MSC-2","MSC-3","Peri-1","Peri-2","Peri-3",
              "Endothelial","osteoblasts","osteoclasts"]
PCT=0.10   # expressed threshold

# ---------- priors ----------
print("loading NicheNet priors...", flush=True)
ltm=pd.read_csv(f"{NN}/ligand_target_matrix.tsv.gz", sep="\t", index_col=0)  # rows=targets, cols=ligands
lr=pd.read_csv(f"{NN}/lr_network.tsv.gz", sep="\t")
lr.columns=[c.lower() for c in lr.columns]
lig_col="from" if "from" in lr.columns else lr.columns[0]
rec_col="to"   if "to"   in lr.columns else lr.columns[1]
print("ltm:", ltm.shape, "| lr_network:", lr.shape, "cols", list(lr.columns), flush=True)

# ---------- receiver (TAM) ----------
print("loading myeloid (receiver)...", flush=True)
mye=sc.read_h5ad(f"{PROJ}/Cleaned_Data/myeloid_FINAL_labels.h5ad")
mye.X=mye.layers["counts"].copy()
tam=mye[mye.obs["final_label"]=="TAM"].copy()
def pct_expressed(ad):
    X=ad.X; X=X.tocsc() if sp.issparse(X) else sp.csc_matrix(X)
    return np.asarray((X>0).mean(axis=0)).ravel()
tam_pct=pd.Series(pct_expressed(tam), index=tam.var_names)
bg_receiver=set(tam_pct[tam_pct>=PCT].index) & set(ltm.index)
print("receiver expressed genes (bg):", len(bg_receiver), flush=True)

# geneset of interest: TAM vs Mono1 up (log_norm)
mye.X=mye.layers["log_norm"].copy()
sub=mye[mye.obs["final_label"].isin(["TAM","Mono1"])].copy()
sc.tl.rank_genes_groups(sub, "final_label", groups=["TAM"], reference="Mono1", method="wilcoxon")
de=sc.get.rank_genes_groups_df(sub, group="TAM")
geneset=de[(de.pvals_adj<0.05)&(de.logfoldchanges>1)]["names"].tolist()
geneset=[g for g in geneset if g in bg_receiver]
pd.DataFrame({"gene":geneset}).to_csv(f"{OUT}/receiver_geneset_TAMvsMono1.csv", index=False)
print("geneset of interest (in bg):", len(geneset), "| LA-TAM in geneset:",
      [g for g in LATAM if g in geneset], flush=True)

# ---------- senders ----------
print("loading integrated (senders)...", flush=True)
ca=pd.read_csv(f"{PROJ}/Cleaned_Data/cell-annotations.csv").set_index("barcode")["celltype"]
integ=sc.read_h5ad(f"{PROJ}/Cleaned_Data/integrated.h5ad")
integ.X=integ.layers["counts"].copy()
lab=ca.reindex(integ.obs_names)
integ.obs["celltype"]=lab.values
send=integ[integ.obs["celltype"].isin(SENDER_TYPES)].copy()
print("sender cells:", send.n_obs, "| per type:\n", send.obs["celltype"].value_counts().to_string(), flush=True)
send_pct=pd.Series(pct_expressed(send), index=send.var_names)
expressed_ligands=set(send_pct[send_pct>=PCT].index)
# per-sender-type expressed ligands (for reporting source)
type_expr={}
for t in SENDER_TYPES:
    st=send[send.obs["celltype"]==t]
    if st.n_obs>=10:
        p=pd.Series(pct_expressed(st), index=st.var_names)
        type_expr[t]=set(p[p>=PCT].index)

# ---------- candidate ligands ----------
lr_expr=lr[(lr[lig_col].isin(expressed_ligands)) & (lr[rec_col].isin(bg_receiver))]
potential_ligands=sorted(set(lr_expr[lig_col]) & set(ltm.columns))
print("potential ligands (expr in sender, receptor in receiver):", len(potential_ligands), flush=True)

# ---------- ligand activity ----------
bg=sorted(bg_receiver)
response=np.array([1.0 if g in set(geneset) else 0.0 for g in bg])
sub_ltm=ltm.loc[bg, potential_ligands]
rows=[]
for L in potential_ligands:
    v=sub_ltm[L].values
    if np.std(v)==0: continue
    r=pearsonr(v, response)[0]
    aupr=average_precision_score(response, v)
    rows.append({"ligand":L, "pearson":r, "aupr":aupr,
                 "sender_types":";".join([t for t in SENDER_TYPES if L in type_expr.get(t,set())])})
act=pd.DataFrame(rows).sort_values("pearson", ascending=False).reset_index(drop=True)
act["rank"]=np.arange(1,len(act)+1)
act.to_csv(f"{TAB}/nichenet_ligand_activity.csv", index=False)
print("top15 ligands by pearson:\n", act.head(15).to_string(), flush=True)

# ---------- ligand -> target links (LA-TAM + NR4A2 focus) ----------
TOPL=act.head(30)["ligand"].tolist()
focus_targets=[g for g in (LATAM+["NR4A2","NR4A3"]) if g in ltm.index]
lt=ltm.loc[focus_targets, TOPL]  # regulatory potential of top ligands -> focus targets
lt.to_csv(f"{TAB}/nichenet_ligand_target_latam.csv")
# long form: ligand-target links limited to geneset targets
links=[]
for L in TOPL:
    col=ltm[L]
    for g in focus_targets:
        links.append({"ligand":L,"target":g,"reg_potential":float(col.get(g,np.nan)),
                      "target_in_geneset":g in set(geneset)})
pd.DataFrame(links).to_csv(f"{TAB}/nichenet_ligand_target_links_long.csv", index=False)

# ---------- ligand -> receptor for top ligands ----------
lr_top=lr_expr[lr_expr[lig_col].isin(TOPL)][[lig_col,rec_col]].drop_duplicates()
lr_top.columns=["ligand","receptor"]
lr_top.to_csv(f"{TAB}/nichenet_ligand_receptor.csv", index=False)

# ---------- figures ----------
plt.figure(figsize=(5,6))
top=act.head(20)[::-1]
plt.barh(top["ligand"], top["pearson"], color="#3b6ea5")
plt.xlabel("Pearson (ligand activity)"); plt.title("Top predicted niche ligands -> LA-TAM program")
plt.tight_layout(); plt.savefig(f"{FIG}/nichenet_ligand_activity.png", dpi=150); plt.close()

plt.figure(figsize=(max(6,0.35*len(TOPL)),0.5*len(focus_targets)+2))
sns.heatmap(lt, cmap="viridis", cbar_kws={"label":"regulatory potential"})
plt.title("Top ligands -> LA-TAM / NR4A targets"); plt.xlabel("ligand"); plt.ylabel("target")
plt.tight_layout(); plt.savefig(f"{FIG}/nichenet_ligand_target_heatmap.png", dpi=150); plt.close()

# NR4A2 downstream summary
nr4a2_rank=act.index[act["ligand"]=="NR4A2"].tolist()
top_for_nr4a2=lt.loc["NR4A2"].sort_values(ascending=False).head(10) if "NR4A2" in lt.index else None
if top_for_nr4a2 is not None:
    top_for_nr4a2.to_csv(f"{TAB}/nichenet_top_ligands_for_NR4A2.csv", header=["reg_potential"])
    print("top ligands predicted to induce NR4A2:\n", top_for_nr4a2.to_string(), flush=True)
print("TASK2 NICHENET DONE", flush=True)
