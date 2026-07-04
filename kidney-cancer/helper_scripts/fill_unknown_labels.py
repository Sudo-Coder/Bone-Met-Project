#!/usr/bin/env python
"""
Fill 'Unknown' cell-type labels by the established label-transfer method:
DISTANCE-WEIGHTED k-NEAREST-NEIGHBOR MAJORITY VOTE in the scVI integration latent
space (X_scVI). This is the standard approach underlying Seurat TransferData /
scanpy ingest / Azimuth, and is more accurate/robust than 1-NN.

- Reference = genuinely annotated cells (from 'cell_type_orig' if present, else the
  current 'cell type'); Query = the originally-'Unknown' cells.
- Writes the predictions back into the ORIGINAL 'cell type' column with NO gaps.
- Preserves the true originals in 'cell_type_orig'; records per-cell confidence
  ('label_transfer_confidence') and mean k-NN distance ('nn_transfer_dist').
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors

PATH = "Cleaned_Data/integrated.h5ad"
LABEL_COL = "cell type"
EMBED = "X_scVI"
UNKNOWN = "Unknown"
K = 15

print(f"[load] {PATH}", flush=True)
a = ad.read_h5ad(PATH)

# True original labels: prefer the preserved 'cell_type_orig' so a prior fill can't
# contaminate the reference set; fall back to the live column on a fresh file.
src_col = "cell_type_orig" if "cell_type_orig" in a.obs.columns else LABEL_COL
orig = a.obs[src_col].astype(str).to_numpy()
Z = np.asarray(a.obsm[EMBED], dtype=np.float32)

query = orig == UNKNOWN          # originally-Unknown cells to fill
ref = ~query                     # genuinely-annotated reference cells
print(f"[info] reference(labeled)={ref.sum()}  query(Unknown)={query.sum()}  total={len(orig)}  embed={Z.shape}", flush=True)
assert query.sum() > 0, "no Unknown cells to fill"

# Distance-weighted k-NN majority vote (the established label-transfer classifier)
clf = KNeighborsClassifier(n_neighbors=K, weights="distance", metric="euclidean", n_jobs=-1)
clf.fit(Z[ref], orig[ref])
pred = clf.predict(Z[query])
conf = clf.predict_proba(Z[query]).max(axis=1)          # vote share of the winning label

# mean distance to the k reference neighbors (transfer-quality audit)
nn = NearestNeighbors(n_neighbors=K, n_jobs=-1).fit(Z[ref])
dists, _ = nn.kneighbors(Z[query])
mean_dist = dists.mean(axis=1)

# ---- write results back into the ORIGINAL column, no gaps ----
new = orig.copy()
new[query] = pred                                        # only Unknown rows change
a.obs["cell_type_orig"] = pd.Categorical(orig)           # preserve/refresh true originals
a.obs[LABEL_COL] = pd.Categorical(new)                   # original column, now gap-free

conf_full = np.ones(len(orig), dtype=np.float32)         # annotated cells = 1.0
conf_full[query] = conf.astype(np.float32)
dist_full = np.zeros(len(orig), dtype=np.float32)
dist_full[query] = mean_dist.astype(np.float32)
a.obs["label_transfer_confidence"] = conf_full
a.obs["nn_transfer_dist"] = dist_full

# ---- hard gap check: every cell has a real label ----
final = a.obs[LABEL_COL].astype(str).to_numpy()
n_unknown = int((final == UNKNOWN).sum())
n_bad = int(pd.isna(a.obs[LABEL_COL]).sum() + (final == "").sum() + (final == "nan").sum())
print(f"[assigned] {query.sum()} filled by k={K} distance-weighted vote", flush=True)
print(pd.Series(pred).value_counts().to_string(), flush=True)
print(f"[confidence] transferred: min={conf.min():.2f} median={np.median(conf):.2f} ; low(<0.5)={int((conf<0.5).sum())}", flush=True)
print(f"[distance]  transferred mean-kNN: median={np.median(mean_dist):.3f} max={mean_dist.max():.3f}", flush=True)
print(f"[GAP CHECK] Unknown_left={n_unknown}  NaN/empty={n_bad}  (both must be 0)", flush=True)
assert n_unknown == 0 and n_bad == 0, "GAPS REMAIN — not saving"

tmp = PATH + ".tmp"
print(f"[save] -> {tmp} then atomic replace", flush=True)
a.write_h5ad(tmp)
os.replace(tmp, PATH)
print("[done] '{}' column filled with no gaps and saved.".format(LABEL_COL), flush=True)
