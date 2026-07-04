#!/usr/bin/env python
"""
Build kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad the SAME way the prostate
myeloid object was made (prostate cellrank c39), from the ATF3-fixed integrated.h5ad.

Recipe:
  1. subset integrated -> myeloid cells = barcodes in myeloid-annotations.csv
  2. myeloid.raw = myeloid                 (full-gene, log-normalized)
  3. HVG (default seurat/dispersion) but FORCE highly_variable=True for ATF3 + AP-1
     so ATF3 survives the cut; then subset to HVGs
  4. scale(max_value=10) -> PCA -> neighbors(n=50, n_pcs=50) -> leiden(res=0.5) -> umap
  5. final_label from myeloid-annotations.csv by barcode
  6. keep condition + Sample + counts + scvi_normalized + leiden (CellOracle needs raw
     counts + a leiden column; CellRank uses the graph / X_scVI)
Output feeds CellRank, then CellOracle.
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc

INTEG = "Cleaned_Data/integrated.h5ad"
ANN   = "Cleaned_Data/myeloid-annotations.csv"
OUT   = "Cleaned_Data/myeloid_FINAL_labels.h5ad"
AP1   = ["ATF3", "JUN", "JUNB", "JUND", "FOS", "FOSB", "EGR1", "NR4A1", "NR4A2", "DUSP1"]

sc.settings.verbosity = 1

print(f"[load] {INTEG}", flush=True)
adata = sc.read_h5ad(INTEG)
print(f"[load] integrated shape {adata.shape}; .X max={float(adata.X[:1000].max()):.2f} (log-norm)", flush=True)

# 1) subset to myeloid using the CORRECTED integrated 'cell type' labels
#    (membership comes from the improved labels; fine labels come from the csv in step 5)
MYELOID_CT = ["Macrophage", "Mono-1", "Mono-2", "Mono-3", "mDC", "Monocyte pro", "osteoclasts"]
ann = pd.read_csv(ANN, index_col="barcode")["celltype"].astype(str)
mask = adata.obs["cell type"].astype(str).isin(MYELOID_CT).to_numpy()
myeloid = adata[mask].copy()
del adata
print(f"[subset] myeloid-by-cell-type: {myeloid.n_obs} cells", flush=True)
print(myeloid.obs["cell type"].astype(str).value_counts().to_string(), flush=True)

# 2) preserve full-gene log-normalized matrix in .raw
myeloid.raw = myeloid

# 3) HVG (default dispersion) + FORCE ATF3/AP-1, then subset
sc.pp.highly_variable_genes(myeloid)                       # flavor='seurat', default cutoffs; flags only
n_auto = int(myeloid.var["highly_variable"].sum())
present = [g for g in AP1 if g in myeloid.var_names]
missing = [g for g in AP1 if g not in myeloid.var_names]
already = [g for g in present if bool(myeloid.var.loc[g, "highly_variable"])]
myeloid.var.loc[present, "highly_variable"] = True         # force AP-1 module in
myeloid = myeloid[:, myeloid.var["highly_variable"]].copy()
print(f"[hvg] auto-HVG={n_auto}; forced AP-1 kept={present} (already-HVG={already}); missing={missing}", flush=True)
print(f"[hvg] final HVG gene count: {myeloid.n_vars}", flush=True)
assert "ATF3" in myeloid.var_names, "ATF3 did not survive HVG subset!"

# 4) recompute on HVGs (prostate cellrank c39): scale -> PCA -> neighbors -> leiden -> umap
sc.pp.scale(myeloid, max_value=10)
sc.pp.pca(myeloid, n_comps=50)
sc.pp.neighbors(myeloid, n_neighbors=50, n_pcs=50)
sc.tl.leiden(myeloid, resolution=0.5, flavor="leidenalg")
sc.tl.umap(myeloid)
print(f"[recompute] leiden clusters: {myeloid.obs['leiden'].nunique()}; obsm={list(myeloid.obsm.keys())}", flush=True)

# 5) impose fine labels from myeloid-annotations.csv by barcode; gaps -> 'Unknown' (kept, per request)
lab = ann.reindex(myeloid.obs_names)
n_unknown = int(lab.isna().sum())
myeloid.obs["final_label"] = pd.Categorical(lab.fillna("Unknown").astype(str).values)
print(f"[final_label] Unknown (myeloid-by-cell-type but not in csv): {n_unknown}", flush=True)
if n_unknown:
    print("[final_label] cell-type composition of the Unknown cells:\n"
          + myeloid.obs.loc[lab.isna().to_numpy(), "cell type"].astype(str).value_counts().to_string(), flush=True)

# 6) reports + sanity for CellOracle/CellRank
print(f"[shape] {myeloid.shape}", flush=True)
print(f"[ATF3 present] {'ATF3' in myeloid.var_names}", flush=True)
print("[final_label counts]\n" + myeloid.obs["final_label"].value_counts().to_string(), flush=True)
print(f"[layers] {list(myeloid.layers.keys())}  (counts = raw for CellOracle)", flush=True)
print(f"[obsm] {list(myeloid.obsm.keys())}", flush=True)
keep_ok = all(c in myeloid.obs.columns for c in ["condition", "Sample", "leiden", "final_label"])
print(f"[obs has condition/Sample/leiden/final_label] {keep_ok}", flush=True)
assert keep_ok and "counts" in myeloid.layers, "missing a required field for downstream tools"

# strip None-valued uns entries: anndata 0.12 (this env) writes them with a 'null'
# encoding that celloracle_env's older anndata 0.10.8 cannot read (e.g. uns['log1p']['base']).
import collections.abc as _abc
def _clean_uns(d):
    for k in list(d.keys()):
        v = d[k]
        if v is None:
            del d[k]
        elif isinstance(v, _abc.MutableMapping):
            _clean_uns(v)
            if len(v) == 0:
                del d[k]
_clean_uns(myeloid.uns)
print(f"[uns] cleaned; keys: {list(myeloid.uns.keys())}", flush=True)

# atomic write
tmp = OUT + ".tmp"
print(f"[save] -> {tmp} then atomic replace", flush=True)
myeloid.write_h5ad(tmp)
os.replace(tmp, OUT)
print(f"[done] wrote {OUT}", flush=True)
