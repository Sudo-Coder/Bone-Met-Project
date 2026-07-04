#!/usr/bin/env python
"""10_label_transfer_withstromal.py — Phase 2-prep: label the RCC with-stromal niche object so CellChat
has sender/receiver identities, and validate that BENIGN samples contain usable stromal/endothelial/
fibroblast senders (the benign-vs-tumor CellChat conditional).

- 81,299/96,829 cells barcode-match the labeled `integrated.h5ad` -> direct 'cell type' transfer.
- the ~15,530 unlabeled cells are the 9 Benign-stroma samples -> kNN-classify in X_scVI (reference
  labeled cells already contain stromal types MSC/Peri/Endo/Osteoclast from tumor samples).
- marker-validate each transferred label; tabulate benign vs tumor sender availability.

Run: envs/rcc_reinterp_venv/bin/python. Seed 0.
"""
import os, warnings, json
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py, anndata as ad
import scipy.sparse as sp
from sklearn.neighbors import KNeighborsClassifier

np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables"); os.makedirs(TAB,exist_ok=True)
WS=os.path.join(ROOT,"kidney-cancer/Cleaned_Data/integrated-with-stromal.h5ad")

MARK={"Tumor_ccRCC":["CA9","PAX8","NDUFA4L2","VEGFA","EPCAM"],
      "Endothelial":["PECAM1","VWF","CLDN5","FLT1"],
      "MSC_Fibroblast":["LUM","DCN","COL1A1","COL1A2","PDGFRA"],
      "Pericyte":["RGS5","ACTA2","PDGFRB","NOTCH3"],
      "Osteoclast":["CTSK","ACP5","MMP9"],
      "Myeloid_TAM":["CD68","LYZ","C1QA","CD14"],
      "T_NK":["CD3D","CD3E","NKG7"],"B":["CD79A","MS4A1"],"Erythroid":["HBB","HBA1"]}

# ---- read with-stromal X (CSR), var, obs, X_scVI via h5py (avoid raw materialization) ----
f=h5py.File(WS,"r")
def idx(g): k=g.attrs.get("_index","_index"); return k.decode() if isinstance(k,bytes) else k
obs=f["obs"]; names=np.array([x.decode() for x in obs[idx(obs)][:]])
var=f["var"]; vnames=np.array([x.decode() for x in var[idx(var)][:]])
sg=obs["Sample"]; scat=np.array([c.decode() for c in sg["categories"][:]]); samp=scat[sg["codes"][:]]
Xg=f["X"]; X=sp.csr_matrix((Xg["data"][:],Xg["indices"][:],Xg["indptr"][:]),shape=tuple(Xg.attrs["shape"]))
scvi=f["obsm"]["X_scVI"][:]
f.close()
gi={g:i for i,g in enumerate(vnames)}

# ---- direct transfer from labeled integrated ----
lab=ad.read_h5ad(os.path.join(ROOT,"kidney-cancer/Cleaned_Data/integrated.h5ad"),backed="r")
lab_ct=lab.obs["cell type"].astype(str); lab_cond=lab.obs["condition"].astype(str)
ct=pd.Series(index=names,dtype=object)
ct_map=lab_ct.reindex(names); ct[:]=ct_map.values
labeled_mask=ct.notna().values
print(f"direct-transferred: {labeled_mask.sum()} / {len(names)}")

# ---- kNN classify the unlabeled (benign-stroma) in X_scVI ----
knn=KNeighborsClassifier(n_neighbors=30,weights="distance")
knn.fit(scvi[labeled_mask], ct[labeled_mask].values)
proba=knn.predict_proba(scvi[~labeled_mask]); classes=knn.classes_
pred=classes[proba.argmax(1)]; conf=proba.max(1)
ct.values[~labeled_mask]=pred
transfer_conf=np.ones(len(names)); transfer_conf[~labeled_mask]=conf
cell_type=ct.astype(str).values

# ---- coarse compartment for CellChat ----
def coarse(c):
    c=str(c)
    if c in ("Tumor",): return "Tumor"
    if c.startswith("MSC"): return "MSC"
    if c.startswith("Peri"): return "Pericyte"
    if c=="Endothelial": return "Endothelial"
    if c in ("Osteoclasts","osteoclasts"): return "Osteoclast"
    if c in ("osteoblasts","Osteoblasts"): return "Osteoblast"
    if c=="TAM": return "TAM"
    if c=="TIM": return "TIM"
    if c in ("Mono1","Mono2","Mono3","Monocyte Pro","mDC","PDC"): return "Myeloid_other"
    if c in ("Thelper","CD4 Naive","CD8 Naive","CTL-1","CTL-2","CTL-3","NKT","NK1","NK2","Treg","Proliferating T"): return "T_NK"
    if "B" in c and c!="Erythroid": return "B"
    if c=="Erythroid": return "Erythroid"
    if c=="Progenitors": return "Progenitors"
    return "other"
comp=np.array([coarse(c) for c in cell_type])
# condition from Sample
def cond_of(s):
    s=str(s)
    if "Benign-stroma" in s: return "Benign_stroma"
    if "Benign-immune" in s: return "Benign_immune"
    if "Tumor" in s: return "Tumor"
    if "Involve" in s: return "Involved"
    if "Noninvolved" in s: return "Distal"
    return "other"
condition=np.array([cond_of(s) for s in samp])

# ---- marker validation: mean log-norm expr per coarse compartment ----
def colmean(gene, mask):
    if gene not in gi: return np.nan
    v=X[:,gi[gene]]; v=v.toarray().ravel()
    return float(v[mask].mean())
rows=[]
for panel,genes in MARK.items():
    for cp in ["Tumor","Endothelial","MSC","Pericyte","Osteoclast","TAM","TIM","T_NK","B","Erythroid"]:
        m=comp==cp
        if m.sum()==0: continue
        vals=[colmean(g,m) for g in genes if g in gi]
        rows.append(dict(marker_panel=panel,compartment=cp,n=int(m.sum()),
                         mean_expr=float(np.nanmean(vals)) if vals else np.nan))
mv=pd.DataFrame(rows); mv.to_csv(os.path.join(TAB,"withstromal_marker_validation.csv"),index=False)

# ---- benign vs tumor sender availability ----
senders=["Tumor","MSC","Pericyte","Endothelial","Osteoclast"]
avail=pd.crosstab(comp,condition)
avail.to_csv(os.path.join(TAB,"withstromal_sender_availability.csv"))

# ---- save labeled metadata for CellChat ----
meta=pd.DataFrame({"cell":names,"Sample":samp,"condition":condition,"cell_type":cell_type,
                   "compartment":comp,"transfer_conf":transfer_conf})
meta.to_csv(os.path.join(TAB,"withstromal_labels.csv"),index=False)

print("\n=== SENDER availability (coarse compartment x condition) ===")
print(avail.reindex(senders+["TAM","TIM","T_NK","B"]).fillna(0).astype(int).to_string())
print("\n=== benign-stroma sender counts (the CellChat benign conditional) ===")
bs=condition=="Benign_stroma"
for cp in senders: print(f"  {cp}: {int((comp[bs]==cp).sum())}")
print("\n=== marker validation (mean log-norm; diagonal should be high) ===")
piv=mv.pivot_table(index="compartment",columns="marker_panel",values="mean_expr")
cols=["Tumor_ccRCC","Endothelial","MSC_Fibroblast","Pericyte","Osteoclast","Myeloid_TAM","T_NK","B","Erythroid"]
print(piv.reindex(columns=[c for c in cols if c in piv.columns]).round(2).to_string())
print("\nlow-confidence transfers (conf<0.5):", int((transfer_conf<0.5).sum()),
      "of", int((~labeled_mask).sum()),"kNN-labeled")
print("wrote withstromal_labels.csv, withstromal_marker_validation.csv, withstromal_sender_availability.csv")
