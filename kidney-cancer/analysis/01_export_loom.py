#!/usr/bin/env python
"""Task 1 step 0 — export raw counts to a loom for pySCENIC, plus an obs metadata table.

Reads (READ-ONLY) the existing myeloid object; writes only into outputs/.
Usage: envs/mechanism_env/bin/python analysis/01_export_loom.py
"""
import os, sys
import numpy as np
import scanpy as sc
import loompy as lp
import pandas as pd
import scipy.sparse as sp

ROOT = "/autofs/projects-t3/hussain/scProj"
PROJ = ROOT + "/kidney-cancer"   # project data/outputs (relocated 2026-07-02)
H5AD = f"{PROJ}/Cleaned_Data/myeloid_FINAL_labels.h5ad"
OUT  = f"{PROJ}/outputs/scenic"
os.makedirs(OUT, exist_ok=True)

print("loading", H5AD, flush=True)
adata = sc.read_h5ad(H5AD)
print("full:", adata.shape, flush=True)

# Use RAW COUNTS (GRNBoost2 expects counts / non-negative expression).
assert "counts" in adata.layers, "counts layer missing"
adata.X = adata.layers["counts"].copy()

# Gene filter: keep genes detected in >= 10 cells (denoise, keep GRN tractable).
# Guarantee genes-of-interest are retained regardless of the filter.
GOI = ["NR4A2","NR4A3","NFKB1","RELA","HIF1A","NFE2L2",
       "TREM2","APOE","APOC1","C1QA","C1QB","C1QC","CD9","FABP5","LIPA","GPNMB"]
X = adata.X
if not sp.issparse(X):
    X = sp.csr_matrix(X)
n_cells_per_gene = np.asarray((X > 0).sum(axis=0)).ravel()
keep = n_cells_per_gene >= 10
for g in GOI:
    if g in adata.var_names:
        keep[adata.var_names.get_loc(g)] = True
adata = adata[:, keep].copy()
print("after gene filter (>=10 cells):", adata.shape, flush=True)
missing = [g for g in GOI if g not in adata.var_names]
print("GOI missing after filter:", missing, flush=True)

# Write loom (cells x genes). pySCENIC expects row=genes attr 'Gene', col=cells attr 'CellID'.
Xc = adata.X
if sp.issparse(Xc):
    Xc = Xc.toarray()
Xc = Xc.astype(np.float32)
row_attrs = {"Gene": np.array(adata.var_names)}
col_attrs = {"CellID": np.array(adata.obs_names)}
loom_path = f"{OUT}/myeloid_counts.loom"
lp.create(loom_path, Xc.T, row_attrs, col_attrs)  # transpose -> genes x cells
print("wrote", loom_path, flush=True)

# Save obs metadata (final_label, condition, sample, GOI expression from log_norm) for downstream joins.
meta = adata.obs[["final_label","condition","Sample","leiden","TAM_score"]].copy()
# also attach log-normalized GOI expression per cell for correlation tests
if "log_norm" in adata.layers:
    ln = adata.layers["log_norm"]
    ln = ln.toarray() if sp.issparse(ln) else np.asarray(ln)
    for g in GOI:
        if g in adata.var_names:
            meta[f"expr_{g}"] = ln[:, adata.var_names.get_loc(g)]
meta.to_parquet(f"{OUT}/cell_metadata.parquet")
print("wrote cell_metadata.parquet with", meta.shape, "cols:", list(meta.columns), flush=True)

# Save the filtered gene universe (background for hypergeometric test in task 1).
pd.Series(list(adata.var_names)).to_csv(f"{OUT}/gene_universe.csv", index=False, header=["gene"])
print("wrote gene_universe.csv:", adata.n_vars, "genes", flush=True)
print("EXPORT DONE", flush=True)
