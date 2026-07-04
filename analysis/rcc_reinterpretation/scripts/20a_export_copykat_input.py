#!/usr/bin/env python
"""20a_export_copykat_input.py — Phase 1.5: export raw counts for CopyKAT malignant-sender validation.
Subset of the labeled RCC full-niche (integrated.h5ad, Tumor condition): Tumor (putative malignant
senders) + stromal senders (MSC/Pericyte/Endothelial) + myeloid (TAM/TIM) + a normal reference (T/NK).
CopyKAT should call Tumor cells aneuploid and stroma/immune diploid -> justifies the malignant sender set.
Run: envs/rcc_reinterp_venv/bin/python. Seed 0.
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, scipy.io as sio
np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables")
OUT=os.path.join(TAB,"copykat_input"); os.makedirs(OUT,exist_ok=True)
a=ad.read_h5ad(os.path.join(ROOT,"kidney-cancer/Cleaned_Data/integrated.h5ad"))
ct=a.obs["cell type"].astype(str); cond=a.obs["condition"].astype(str)
def coarse(c):
    if c=="Tumor": return "Tumor"
    if c.startswith("MSC"): return "MSC"
    if c.startswith("Peri"): return "Pericyte"
    if c=="Endothelial": return "Endothelial"
    if c in ("TAM","TIM"): return "Myeloid_TAM"
    if c in ("Thelper","CD4 Naive","CD8 Naive","CTL-1","CTL-2","CTL-3","NKT","NK1","NK2","Treg"): return "T_NK"
    return "other"
comp=ct.map(coarse)
tum=cond.eq("Tumor").values
rng=np.random.default_rng(0)
pick=[]
for grp,n in [("Tumor",99999),("MSC",99999),("Pericyte",99999),("Endothelial",99999),("Myeloid_TAM",800),("T_NK",1500)]:
    idx=np.where((comp.values==grp)&tum)[0]
    if len(idx)>n: idx=rng.choice(idx,n,replace=False)
    pick.append(idx)
pick=np.concatenate(pick)
sub=a[pick]
counts=sub.layers["counts"]; counts=counts if sp.issparse(counts) else sp.csr_matrix(counts)
sio.mmwrite(os.path.join(OUT,"counts.mtx"), counts.T.tocsc())   # genes x cells
pd.Series(sub.var_names).to_csv(os.path.join(OUT,"genes.txt"),index=False,header=False)
pd.Series(sub.obs_names).to_csv(os.path.join(OUT,"cells.txt"),index=False,header=False)
meta=pd.DataFrame({"cell":sub.obs_names,"compartment":comp.values[pick],"cell_type":ct.values[pick]})
meta.to_csv(os.path.join(OUT,"meta.csv"),index=False)
# known-normal reference = T_NK barcodes
meta.loc[meta.compartment=="T_NK","cell"].to_csv(os.path.join(OUT,"normal_cells.txt"),index=False,header=False)
print("exported", sub.n_obs, "cells x", sub.n_vars, "genes")
print(meta["compartment"].value_counts().to_string())
