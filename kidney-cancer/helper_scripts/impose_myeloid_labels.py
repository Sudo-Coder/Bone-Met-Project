#!/usr/bin/env python
"""
Impose the refined myeloid final_label onto integrated.h5ad.
- drop the 40 myeloid-by-cell-type cells that were removed from the myeloid file (cluster 15)
- back up integrated 'cell type' -> 'cell_type_coarse'
- overwrite 'cell type': myeloid cells (by barcode) get their fine final_label,
  everyone else keeps their existing label
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import collections.abc as abc

INTEG = "Cleaned_Data/integrated.h5ad"
MYE   = "Cleaned_Data/myeloid_FINAL_labels.h5ad"
MYELOID_CT = ["Macrophage", "Mono-1", "Mono-2", "Mono-3", "mDC", "Monocyte pro", "osteoclasts"]

# barcode -> fine myeloid label
mye = sc.read_h5ad(MYE, backed="r")
fine = mye.obs["final_label"].astype(str)          # indexed by barcode
mye_bc = set(fine.index)
print(f"[myeloid] {len(mye_bc)} cells; labels: {dict(fine.value_counts())}", flush=True)

print(f"[load] {INTEG}", flush=True)
a = sc.read_h5ad(INTEG)
before = a.n_obs

# 1) drop myeloid-by-cell-type cells NOT in the myeloid file (the removed cluster-15 cells)
is_mye_ct = a.obs["cell type"].astype(str).isin(MYELOID_CT).to_numpy()
in_file   = a.obs_names.isin(mye_bc)
drop_mask = is_mye_ct & (~in_file)
print(f"[drop] {int(drop_mask.sum())} cells (myeloid-by-cell-type, not in myeloid file)", flush=True)
a = a[~drop_mask].copy()

# 2) back up coarse labels
a.obs["cell_type_coarse"] = a.obs["cell type"].astype(str).values

# 3) overwrite 'cell type' with the fine myeloid label where available
fine_full = pd.Series(a.obs_names, index=a.obs_names).map(fine)   # NaN for non-myeloid
ct = a.obs["cell type"].astype(str).copy()
mask = fine_full.notna().to_numpy()
ct[mask] = fine_full[mask].astype(str).values
a.obs["cell type"] = pd.Categorical(ct)

# reports / sanity
print(f"[result] {before} -> {a.n_obs} cells; imposed fine label on {int(mask.sum())} myeloid cells", flush=True)
left = int(a.obs["cell type"].astype(str).isin(MYELOID_CT).sum())
print(f"[check] coarse myeloid labels remaining in 'cell type': {left} (should be 0)", flush=True)
print(f"[check] NaN/Unknown in 'cell type': "
      f"{int(pd.isna(a.obs['cell type']).sum())} / {int((a.obs['cell type'].astype(str)=='Unknown').sum())}", flush=True)
print("[cell type value_counts]\n" + a.obs["cell type"].value_counts().to_string(), flush=True)

# clean None-valued uns (anndata 0.10.8 / celloracle compat) + atomic save
def clean(d):
    for k in list(d.keys()):
        v = d[k]
        if v is None:
            del d[k]
        elif isinstance(v, abc.MutableMapping):
            clean(v)
            if len(v) == 0:
                del d[k]
clean(a.uns)

tmp = INTEG + ".tmp"
print(f"[save] -> {tmp} then atomic replace", flush=True)
a.write_h5ad(tmp)
os.replace(tmp, INTEG)
print(f"[done] wrote {INTEG}", flush=True)
