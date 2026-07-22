#!/usr/bin/env python
# Fair comparison for the 25-gene signature.
#
# 08_prognostic_model.py compares a FITTED LASSO-Cox model against UNFITTED mean-z-score gene
# sets (Obradovic, panTAM) and against random genes with RANDOM weights. A fitted model will
# beat an unfitted score almost by construction, so that comparison does not isolate the
# biology. Here every comparator goes through the identical pipeline -- same train/test split,
# same univariate p<0.05 prefilter, same LASSO-Cox with 10-fold CV -- so the only thing that
# differs is which genes were available.
import os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored
from sklearn.model_selection import train_test_split, KFold
np.random.seed(0)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
MOD = os.path.join(ROOT, "analysis/rcc_reinterpretation/outputs/model")
N_RANDOM = int(os.environ.get("N_RANDOM", "100"))

cand = pd.read_csv(os.path.join(MOD, "candidate_genes.csv"))
kirc = pd.read_csv(os.path.join(ROOT, "resources/tcga/KIRC_HiSeqV2.gz"), sep = "\t", index_col = 0)
kirc.index = [g.upper() for g in kirc.index]
surv = pd.read_csv(os.path.join(ROOT, "resources/tcga/KIRC_survival.txt"), sep = "\t")
surv["bcr"] = surv["sample"].str[:15]

pool = [g for g in cand["gene"] if g in kirc.index]
expr = kirc.loc[pool].T
expr["bcr"] = [s[:15] for s in expr.index]
clin = surv[["bcr", "OS", "OS.time"]].dropna()
clin = clin[clin["OS.time"] > 0]
df = expr.merge(clin, on = "bcr")

# identical split to 08_prognostic_model.py
tr_idx, te_idx = train_test_split(np.arange(len(df)), test_size = 0.3, random_state = 0,
                                  stratify = df["OS"])
train = df.iloc[tr_idx].reset_index(drop = True)
test = df.iloc[te_idx].reset_index(drop = True)
y_tr = Surv.from_arrays(train["OS"].astype(bool), train["OS.time"])
ev_te = test["OS"].astype(bool).values
t_te = test["OS.time"].values

def fit_eval(genes, prefilter = True):
    """same pipeline as the real signature: univariate screen -> LASSO-Cox -> held-out C-index"""
    g = [x for x in genes if x in train.columns]
    if len(g) < 2: return np.nan, 0
    if prefilter:
        keep = []
        for q in g:
            d = train[[q, "OS", "OS.time"]].copy()
            d["z"] = (d[q] - d[q].mean()) / d[q].std()
            try:
                r = CoxPHFitter().fit(d[["z", "OS.time", "OS"]], "OS.time", "OS").summary.loc["z"]
                if r["p"] < 0.05: keep.append(q)
            except Exception:
                pass
        g = keep if len(keep) >= 2 else g
    mu, sd = train[g].mean(), train[g].std()
    xz = ((train[g] - mu) / sd).values
    try:
        path = CoxnetSurvivalAnalysis(l1_ratio = 1.0, alpha_min_ratio = 0.01, max_iter = 100000).fit(xz, y_tr)
    except Exception:
        return np.nan, len(g)
    kf = KFold(n_splits = 10, shuffle = True, random_state = 0)
    cv = []
    for a in path.alphas_:
        fold = []
        for i_tr, i_te in kf.split(xz):
            try:
                m = CoxnetSurvivalAnalysis(l1_ratio = 1.0, alphas = [a], max_iter = 100000).fit(xz[i_tr], y_tr[i_tr])
                fold.append(concordance_index_censored(train["OS"].astype(bool).values[i_te],
                                                       train["OS.time"].values[i_te],
                                                       m.predict(xz[i_te]))[0])
            except Exception:
                fold.append(np.nan)
        cv.append(np.nanmean(fold))
    best = path.alphas_[int(np.nanargmax(cv))]
    m = CoxnetSurvivalAnalysis(l1_ratio = 1.0, alphas = [best], max_iter = 100000).fit(xz, y_tr)
    nz = int((m.coef_ != 0).sum())
    risk = m.predict(((test[g] - mu) / sd).values)
    return float(concordance_index_censored(ev_te, t_te, risk)[0]), nz

rows = []
model = json.load(open(os.path.join(MOD, "risk_model.json")))
mu, sd, co = np.array(model["mean"]), np.array(model["sd"]), np.array(model["coef"])
gg = list(model["genes"])
risk = ((test[gg] - mu) / sd).values @ co
rows.append(dict(comparator = "25-gene signature (frozen)", n_genes_in = len(gg), n_selected = len(gg),
                 heldout_cindex = float(concordance_index_censored(ev_te, t_te, risk)[0])))

for nm, gs in [("Obradovic TREM2 genes", ["TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "GPNMB",
                                          "FOLR2", "SPP1", "CTSD", "CD68"]),
               ("panTAM genes", ["CD68", "CD163", "MRC1", "CSF1R", "LYZ", "AIF1", "FCGR3A"]),
               ("CLEC_LAM8 core", ["C1QA", "C1QB", "C1QC", "APOE", "APOC1", "TREM2", "GPNMB", "MERTK"])]:
    c, nz = fit_eval(gs)
    rows.append(dict(comparator = nm + " (refit, same pipeline)", n_genes_in = len(gs),
                     n_selected = nz, heldout_cindex = c))

rng = np.random.default_rng(0)
rc = []
for i in range(N_RANDOM):
    gs = list(rng.choice(pool, size = len(gg), replace = False))
    c, _ = fit_eval(gs, prefilter = True)
    if c == c: rc.append(c)
    if (i + 1) % 20 == 0: print("  random refits: %d/%d" % (i + 1, N_RANDOM), flush = True)
rc = np.array(rc)
rows.append(dict(comparator = "random %d-gene sets from same pool (refit), mean" % len(gg),
                 n_genes_in = len(gg), n_selected = np.nan, heldout_cindex = float(rc.mean())))
rows.append(dict(comparator = "random sets (refit), 95th percentile",
                 n_genes_in = len(gg), n_selected = np.nan, heldout_cindex = float(np.percentile(rc, 95))))

pd.DataFrame({"cindex": rc}).to_csv(os.path.join(MOD, "fair_benchmark_random_draws.csv"), index = False)
R = pd.DataFrame(rows)
sig_c = R.loc[0, "heldout_cindex"]
R["vs_signature"] = (sig_c - R["heldout_cindex"]).round(3)
R.to_csv(os.path.join(MOD, "fair_benchmark.csv"), index = False)
emp_p = float((rc >= sig_c).mean())
pd.set_option("display.width", 200)
print("\n=== held-out C-index, every comparator through the same fitting pipeline ===")
print(R.round(3).to_string(index = False))
print("\nempirical p (random refits >= signature): %.3f  (n=%d draws)" % (emp_p, len(rc)))
print("random-draw distribution: min %.3f | 25th %.3f | median %.3f | 75th %.3f | max %.3f"
      % (rc.min(), np.percentile(rc, 25), np.median(rc), np.percentile(rc, 75), rc.max()))
for thr in [0.65, 0.70, 0.75]:
    print("  fraction of random draws with C-index > %.2f: %.1f%%" % (thr, 100*(rc > thr).mean()))
