#!/usr/bin/env python
import os, warnings, json
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, decoupler as dc
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
np.random.seed(0)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
TAB = os.path.join(ROOT, "analysis", "rcc_reinterpretation", "outputs", "tables")
os.makedirs(TAB, exist_ok = True)
MIN = 20
MODULES = {
 "CLEC_LAM8": ["C1QA", "C1QB", "C1QC", "APOE", "APOC1", "TREM2", "GPNMB", "MERTK"],
 "RCC_skew_CORE": ["C1QA", "C1QB", "C1QC", "APOE", "APOC1", "TREM2"],
 "complement_C1Q": ["C1QA", "C1QB", "C1QC"],
 "complement_C1Q_C3": ["C1QA", "C1QB", "C1QC", "C3"],
 "APOE_TREM2": ["APOE", "APOC1", "TREM2", "TYROBP"],
 "MERTK_GPNMB_support": ["MERTK", "GPNMB"],
 "SPP1_TAM": ["SPP1", "FN1", "MMP9", "CTSB", "CTSD"],
 "inflammatory_mono": ["IL1B", "CXCL8", "S100A8", "S100A9", "FCN1", "VCAN"],
 "panTAM": ["CD68", "CD163", "MRC1", "CSF1R", "LYZ"],
}
OBJS = {"RCC": dict(path = "kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad", lab = "final_label"),
      "prostate": dict(path = "prostate-cancer/Cleaned_Data/myeloid_FINAL.h5ad", lab = "scanvi_labels_annotation_model_refined")}
def parse_patient(s):
    s = str(s).split("_", 1)[-1].replace(".count", "")
    for kw in ["-Tumor", "-Involve", "-Noninvolved", "-Distal", "-Benign"]: s = s.split(kw)[0]
    return s.replace("RCC-", "")
def norm_cond(c): return {"Involve": "Involved", "Noninvolved": "Distal"}.get(str(c), str(c))

per = []
for tag, cfg in OBJS.items():
    a = ad.read_h5ad(os.path.join(ROOT, cfg["path"]))
    s = a.raw.to_adata() if (a.raw is not None and a.raw.n_vars>=a.n_vars) else a.copy()
    s.obs = a.obs.copy()
    nets = []
    for m, g in MODULES.items():
        gg = [x for x in g if x in s.var_names]
        nets.append(pd.DataFrame({"source": m, "target": gg, "weight": 1.0}))
    net = pd.concat(nets)
    dc.mt.aucell(s, net, tmin = 2, verbose = False)
    sc = s.obsm["score_aucell"]
    df = pd.DataFrame(index = a.obs_names)
    for m in MODULES: df[m] = sc[m].values if m in sc.columns else np.nan
    df["compartment"] = a.obs[cfg["lab"]].astype(str).values
    df["condition"] = a.obs["condition"].astype(str).map(norm_cond).values
    df["patient_id"] = a.obs["Sample"].map(parse_patient).values
    df["cancer_type"] = tag
    df["dataset"] = "rcc_bm" if tag=="RCC" else "prostate_kfoury"
    df["Sample"] = a.obs["Sample"].values
    per.append(df)
pc = pd.concat(per)
pc.to_parquet(os.path.join(TAB, "module_scores_percell.parquet"))

def pseudobulk(df):
    g = df.groupby(["Sample", "cancer_type", "dataset", "patient_id", "condition"], observed = True)
    return g.agg(**{**{m: (m, "mean") for m in MODULES}, "n_cells": ("condition", "size")}).reset_index()
pb_all = pseudobulk(pc)
pb_all = pb_all[pb_all.n_cells>=MIN]
pb_tt = pseudobulk(pc[pc.compartment.isin(["TAM", "TIM"])])
pb_tt = pb_tt[pb_tt.n_cells>=10]
pb_all.to_csv(os.path.join(TAB, "module_pseudobulk_all_myeloid.csv"), index = False)
pb_tt.to_csv(os.path.join(TAB, "module_pseudobulk_TAMTIM.csv"), index = False)

form = "{m} ~ C(cancer_type, Treatment('prostate'))*C(condition, Treatment('Benign'))"
coef = "C(cancer_type, Treatment('prostate'))[T.RCC]:C(condition, Treatment('Benign'))[T.Tumor]"
def interaction(pb, m):
    d = pb[pb.condition.isin(["Tumor", "Benign"])].copy()
    d = d[d[m].notna()]
    if len(d)<6 or d[m].nunique()<3:
        return dict(module = m, unit = "all_myeloid" if pb is pb_all else "TAM+TIM", method = "skip", estimate = np.nan, ci_low = np.nan, ci_high = np.nan, p = np.nan, score_SD = np.nan, est_SD = np.nan, ci_low_SD = np.nan, ci_high_SD = np.nan)
    d["condition"] = pd.Categorical(d["condition"], ["Benign", "Tumor"])
    d["cancer_type"] = pd.Categorical(d["cancer_type"], ["prostate", "RCC"])
    sd = float(d[m].std())
    try:
        mm = smf.mixedlm(form.format(m = m), d, groups = d["patient_id"]).fit(reml = False, method = "lbfgs")
        est, se, p = mm.params[coef], mm.bse[coef], mm.pvalues[coef]
        meth = "mixedlm"
    except Exception:
        mm = smf.ols(form.format(m = m), d).fit(cov_type = "cluster", cov_kwds = {"groups": d["patient_id"]})
        ci = mm.conf_int().loc[coef]
        est, se, p = mm.params[coef], (ci[1]-ci[0])/3.92, mm.pvalues[coef]
        meth = "ols_cluster"
    lo, hi = est-1.96*se, est+1.96*se
    return dict(module = m, unit = "all_myeloid" if pb is pb_all else "TAM+TIM", method = meth,
                estimate = est, ci_low = lo, ci_high = hi, p = p, score_SD = sd,
                est_SD = est/sd, ci_low_SD = lo/sd, ci_high_SD = hi/sd)
rows = []
for m in MODULES:
    rows.append(interaction(pb_all, m))
    rows.append(interaction(pb_tt, m))
res = pd.DataFrame(rows)

prim = res[res.unit=="all_myeloid"].copy()
prim["fdr"] = multipletests(prim["p"], method = "fdr_bh")[1]
res = res.merge(prim[["module", "unit", "fdr"]], on = ["module", "unit"], how = "left")
res.to_csv(os.path.join(TAB, "module_interactions.csv"), index = False)

pd.set_option("display.width", 220, "display.max_colwidth", 22)
print("=== SCORE SCALE ===\nAUCell rank-based enrichment (0..1); NOT z-scored, NOT raw log-norm mean.")
print("Pseudobulk all-myeloid score SD (the estimand) per module:")
for m in MODULES:
    r = res[(res.module==m)&(res.unit=='all_myeloid')].iloc[0]
    print(f"  {m:22s} SD={r.score_SD:.4f}")
print("\n=== CROSS-CANCER INTERACTION (all-myeloid PRIMARY) — raw and SD units ===")
show = res[res.unit=="all_myeloid"][["module", "estimate", "ci_low", "ci_high", "p", "fdr", "est_SD", "ci_low_SD", "ci_high_SD"]]
print(show.round(4).to_string(index = False))
print("\n=== TAM+TIM secondary ===")
print(res[res.unit=="TAM+TIM"][["module", "estimate", "ci_low", "ci_high", "p", "est_SD"]].round(4).to_string(index = False))
