#!/usr/bin/env python
"""Task 4 — recompute a Mono1-rooted pseudotime over the myeloid object and test whether
NR4A2/NR4A3 rise BEFORE the LA-TAM program along the Mono1->TAM lineage.

Approach (reproducible; the CellRank notebook's Palantir pseudotime is never saved to disk):
  - neighbors on X_scVI -> diffmap -> DPT rooted at a Mono1 cell (min TAM_score).
  - CellRank PseudotimeKernel -> GPCCA -> TAM fate probability as lineage weight (best-effort).
  - lineage = monocytes + TAM; smoothed gene trends (lowess) vs pseudotime.
  - per-gene onset (half-max) + peak pseudotime; ordering NR4A vs LA-TAM; cross-correlation lag.
Run in cellrank_env:
  envs/cellrank_env/bin/python analysis/07_pseudotime_trends.py
Writes only outputs/ (pseudotime.parquet consumed by task 3 overlay).
"""
import os, numpy as np, pandas as pd, scanpy as sc
import scipy.sparse as sp
from statsmodels.nonparametric.smoothers_lowess import lowess
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT="/autofs/projects-t3/hussain/scProj"
PROJ = ROOT + "/kidney-cancer"   # project data/outputs (relocated 2026-07-02)
TAB=f"{PROJ}/outputs/tables"; FIG=f"{PROJ}/outputs/figures"; CR=f"{PROJ}/outputs/cellrank"
for d in (TAB,FIG,CR): os.makedirs(d, exist_ok=True)
LATAM=["TREM2","APOE","APOC1","C1QA","C1QB","C1QC","CD9","FABP5"]
NR4A=["NR4A2","NR4A3"]

print("loading myeloid...", flush=True)
ad=sc.read_h5ad(f"{PROJ}/Cleaned_Data/myeloid_FINAL_labels.h5ad")
ad.obs["final_label"]=ad.obs["final_label"].astype("category")   # keep categorical for CellRank .cat

# ---- neighbors on scVI latent, diffmap, DPT rooted at Mono1 ----
sc.pp.neighbors(ad, use_rep="X_scVI", n_neighbors=30)
sc.tl.diffmap(ad)
mono1=np.where(ad.obs["final_label"].values=="Mono1")[0]
root=mono1[np.argmin(ad.obs["TAM_score"].values[mono1])]   # most naive Mono1
ad.uns["iroot"]=int(root)
sc.tl.dpt(ad)
ad.obs["dpt_pseudotime"]=ad.obs["dpt_pseudotime"].replace([np.inf,-np.inf],np.nan)
print("DPT done; root cell:", ad.obs_names[root], flush=True)

# ---- CellRank TAM fate probability (best effort) ----
tam_fate=None
try:
    import cellrank as cr
    pk=cr.kernels.PseudotimeKernel(ad, time_key="dpt_pseudotime").compute_transition_matrix()
    g=cr.estimators.GPCCA(pk); g.compute_schur(n_components=10)
    g.compute_macrostates(n_states=6, cluster_key="final_label")
    g.predict_terminal_states()
    g.compute_fate_probabilities()
    fp=g.fate_probabilities
    names=list(fp.names)
    tam_cols=[i for i,n in enumerate(names) if "TAM" in n]
    if tam_cols:
        tam_fate=pd.Series(np.asarray(fp)[:,tam_cols].sum(1), index=ad.obs_names, name="TAM_fate")
        print("TAM fate prob computed; terminal states:", names, flush=True)
except Exception as e:
    print("CellRank fate step skipped:", type(e).__name__, str(e)[:200], flush=True)

# save pseudotime (+fate) for downstream overlays
pt=ad.obs[["dpt_pseudotime","final_label","condition"]].copy()
if tam_fate is not None: pt["TAM_fate"]=tam_fate
pt.to_csv(f"{CR}/pseudotime.csv")                 # CSV (cellrank_env has no pyarrow)
print("wrote pseudotime.csv", flush=True)

# ---- lineage cells = monocytes + TAM (the Mono1->TAM axis) ----
lin_types=["Mono1","Mono2","Mono3","Monocyte Pro","TAM"]
mask=ad.obs["final_label"].isin(lin_types).values & np.isfinite(ad.obs["dpt_pseudotime"].values)
lin=ad[mask].copy()
ps=lin.obs["dpt_pseudotime"].values
# expression (log_norm)
Xln=lin.layers["log_norm"]; Xln=Xln.toarray() if sp.issparse(Xln) else np.asarray(Xln)
genes=[g for g in NR4A+LATAM if g in lin.var_names]
gi={g:lin.var_names.get_loc(g) for g in genes}
order=np.argsort(ps); psx=ps[order]

def smooth(y):
    z=lowess(y[order], psx, frac=0.3, return_sorted=False)
    return z
def norm01(z):
    z=np.asarray(z); lo,hi=np.nanmin(z),np.nanmax(z)
    return (z-lo)/(hi-lo+1e-12)

trends={g:smooth(Xln[:,gi[g]]) for g in genes}
# onset (first pseudotime reaching 50% of max of smoothed, normalized) and peak
rows=[]
for g in genes:
    z=norm01(trends[g]);
    peak=psx[int(np.nanargmax(z))]
    onset_idx=np.where(z>=0.5)[0]
    onset=psx[onset_idx[0]] if len(onset_idx) else np.nan
    rows.append({"gene":g,"group":"NR4A" if g in NR4A else "LA-TAM",
                 "onset_pt_halfmax":float(onset),"peak_pt":float(peak)})
pk_df=pd.DataFrame(rows).sort_values("onset_pt_halfmax")
pk_df.to_csv(f"{TAB}/genetrend_peaktimes.csv", index=False)
print(pk_df.to_string(), flush=True)

# summary ordering test: mean NR4A onset vs mean LA-TAM onset
nr4a_onset=pk_df[pk_df.group=="NR4A"]["onset_pt_halfmax"].mean()
latam_onset=pk_df[pk_df.group=="LA-TAM"]["onset_pt_halfmax"].mean()
# cross-correlation lag between mean NR4A trend and mean LA-TAM trend
mean_nr4a=norm01(np.mean([norm01(trends[g]) for g in NR4A if g in trends],axis=0))
mean_latam=norm01(np.mean([norm01(trends[g]) for g in LATAM if g in trends],axis=0))
lags=np.arange(-len(psx)+1,len(psx))
xcorr=np.correlate(mean_nr4a-mean_nr4a.mean(), mean_latam-mean_latam.mean(), mode="full")
best_lag=lags[np.argmax(xcorr)]   # >0 => NR4A leads LA-TAM
summary=pd.DataFrame([{"mean_NR4A_onset":nr4a_onset,"mean_LATAM_onset":latam_onset,
   "NR4A_precedes_LATAM": bool(nr4a_onset<latam_onset),
   "xcorr_best_lag_cells_NR4A_leads": int(best_lag)}])
summary.to_csv(f"{TAB}/genetrend_ordering_summary.csv", index=False)
print(summary.to_string(), flush=True)

# ---- figure ----
plt.figure(figsize=(8,5))
for g in NR4A:
    if g in trends: plt.plot(psx, norm01(trends[g]), lw=2.5, ls="--", label=g)
for g in LATAM:
    if g in trends: plt.plot(psx, norm01(trends[g]), lw=1.5, alpha=0.8, label=g)
plt.xlabel("DPT pseudotime (Mono1 -> TAM)"); plt.ylabel("smoothed expression (0-1)")
plt.title("NR4A vs LA-TAM program along the Mono1->TAM lineage")
plt.legend(ncol=2, fontsize=8); plt.tight_layout()
plt.savefig(f"{FIG}/genetrends_nr4a_vs_latam.png", dpi=150); plt.close()
print("TASK4 PSEUDOTIME TRENDS DONE", flush=True)
