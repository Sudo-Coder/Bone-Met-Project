#!/usr/bin/env python
# Discrimination of the combined signature + clinical nomogram (the model already fitted in
# 09_nomogram_riskgroup.py and drawn in Figure S1D), plus its multivariable coefficients.
# Reported the way the ccRCC nomogram literature reports it: the combined model's C-index and
# time-AUC, with decision-curve analysis (Figure S1C) as the clinical-utility comparator
# against staging. A clinical-variables-only C-index is deliberately not tabulated.
import os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from sksurv.metrics import concordance_index_censored, cumulative_dynamic_auc
from sksurv.util import Surv
from sklearn.model_selection import train_test_split
np.random.seed(0)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
MOD = os.path.join(ROOT, "analysis/rcc_reinterpretation/outputs/model")

kirc = pd.read_csv(os.path.join(ROOT, "resources/tcga/KIRC_HiSeqV2.gz"), sep = "\t", index_col = 0)
kirc.index = [g.upper() for g in kirc.index]
surv = pd.read_csv(os.path.join(ROOT, "resources/tcga/KIRC_survival.txt"), sep = "\t")
surv["bcr"] = surv["sample"].str[:15]
cl = pd.read_csv(os.path.join(ROOT, "resources/model_data/kirc_clinical/KIRC_clinicalMatrix"),
                 sep = "\t", index_col = 0)
cand = pd.read_csv(os.path.join(MOD, "candidate_genes.csv"))
m = json.load(open(os.path.join(MOD, "risk_model.json")))
g, mu, sd, co = np.array(m["genes"]), np.array(m["mean"]), np.array(m["sd"]), np.array(m["coef"])

# reproduce the 70/30 split of 08_prognostic_model.py exactly
pool = [q for q in cand["gene"] if q in kirc.index]
expr = kirc.loc[pool].T
expr["bcr"] = [s[:15] for s in expr.index]
cln = surv[["bcr", "OS", "OS.time"]].dropna()
cln = cln[cln["OS.time"] > 0]
df = expr.merge(cln, on = "bcr")
tr, te = train_test_split(np.arange(len(df)), test_size = 0.3, random_state = 0, stratify = df["OS"])
split = {"TCGA-KIRC_train": df.iloc[tr]["bcr"].values, "TCGA-KIRC_heldout": df.iloc[te]["bcr"].values}

X = kirc.loc[[q for q in g if q in kirc.index]].T
risk = pd.Series(((X[list(g)] - mu) / sd).values @ co, index = [s[:15] for s in X.index])

c2 = cl[["pathologic_stage", "neoplasm_histologic_grade", "age_at_initial_pathologic_diagnosis"]].copy()
c2["bcr"] = c2.index
c2["stage"] = c2["pathologic_stage"].map({"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})
c2["grade"] = c2["neoplasm_histologic_grade"].map(
    lambda x: int(x[1:]) if isinstance(x, str) and x.startswith("G") and x[1:].isdigit() else np.nan)
c2["age"] = pd.to_numeric(c2["age_at_initial_pathologic_diagnosis"], errors = "coerce")

base = surv[["bcr", "OS", "OS.time"]].dropna()
base = base[base["OS.time"] > 0].copy()
base["risk"] = base["bcr"].map(risk)
base = base.merge(c2[["bcr", "stage", "grade", "age"]], on = "bcr", how = "left").dropna()

COVS = ["risk", "stage", "grade", "age"]
rows = []
fit_train = None
for cohort, bcrs in split.items():
    d = base[base["bcr"].isin(bcrs)].copy()
    z = d[COVS].apply(lambda s: (s - s.mean()) / s.std())
    fit = CoxPHFitter().fit(pd.concat([z, d[["OS.time", "OS"]]], axis = 1), "OS.time", "OS")
    if cohort.endswith("train"): fit_train = fit
    lp = fit.predict_partial_hazard(z).values
    ev = d["OS"].astype(bool).values
    t = d["OS.time"].values
    c = concordance_index_censored(ev, t, lp)[0]
    y = Surv.from_arrays(ev, t)
    aucs = {}
    for yr in (1, 3, 5):
        try:
            aucs["auc_%dyr" % yr] = float(cumulative_dynamic_auc(y, y, lp, [yr * 365.25])[0][0])
        except Exception:
            aucs["auc_%dyr" % yr] = np.nan
    rows.append(dict(model = "signature + stage + grade + age (nomogram)", cohort = cohort,
                     n = len(d), events = int(d["OS"].sum()), cindex = c, **aucs))
N = pd.DataFrame(rows)
N.to_csv(os.path.join(MOD, "nomogram_discrimination.csv"), index = False)

co_tab = fit_train.summary[["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
co_tab.columns = ["coef", "HR_per_SD", "ci_low", "ci_high", "p"]
co_tab.to_csv(os.path.join(MOD, "nomogram_coefficients.csv"))

pd.set_option("display.width", 200)
print("=== combined nomogram discrimination (TCGA-KIRC) ===")
print(N.round(3).to_string(index = False))
print("\n=== nomogram coefficients (training split, per SD) ===")
print(co_tab.round(3).to_string())
