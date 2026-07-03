#!/usr/bin/env python
"""Task 6 (optional) — TCGA-KIRC external validation.
Score the LA-TAM (CD9+) signature and NR4A2 in TCGA-KIRC bulk RNA-seq (UCSC Xena),
test association with overall survival (Cox + KM by tertile) and correlation LA-TAM vs NR4A2.
Run: envs/mechanism_env/bin/python analysis/08_tcga_kirc.py
Writes only outputs/.
"""
import os, io, sys, numpy as np, pandas as pd, urllib.request
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from scipy.stats import spearmanr, zscore

ROOT="/autofs/projects-t3/hussain/scProj"
PROJ = ROOT + "/kidney-cancer"   # project data/outputs (relocated 2026-07-02)
TAB=f"{PROJ}/outputs/tables"; FIG=f"{PROJ}/outputs/figures"; RES=f"{ROOT}/resources/tcga"
for d in (TAB,FIG,RES): os.makedirs(d, exist_ok=True)
LATAM=["TREM2","APOE","APOC1","C1QA","C1QB","C1QC","CD9","FABP5"]

EXPR_URL="https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/HiSeqV2.gz"          # log2(RSEM+1), gene symbols
SURV_URLS=["https://tcga.xenahubs.net/download/survival/KIRC_survival.txt",
           "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FKIRC_survival.txt.gz"]

def fetch(url, dst):
    if os.path.exists(dst) and os.path.getsize(dst)>0: return dst
    print("downloading", url, flush=True)
    urllib.request.urlretrieve(url, dst); return dst

# expression
expr_path=fetch(EXPR_URL, f"{RES}/KIRC_HiSeqV2.gz")
expr=pd.read_csv(expr_path, sep="\t", index_col=0)   # genes x samples
print("expr:", expr.shape, flush=True)

# survival
surv=None
for u in SURV_URLS:
    try:
        p=fetch(u, f"{RES}/KIRC_survival"+(".gz" if u.endswith(".gz") else ".txt"))
        surv=pd.read_csv(p, sep="\t")
        break
    except Exception as e:
        print("surv url failed:", u, e, flush=True)
if surv is None:
    print("NO SURVIVAL DATA - aborting task6 gracefully", flush=True); sys.exit(0)
print("surv cols:", list(surv.columns)[:10], flush=True)
# standardize survival columns (Xena: sample, OS, OS.time)
scol=[c for c in surv.columns if c.lower() in ("sample","samples","sampleid","_patient")][0]
os_ev=[c for c in surv.columns if c.upper() in ("OS","OS_EVENT","_EVENT")][0]
os_t =[c for c in surv.columns if c.upper() in ("OS.TIME","OS_TIME","_TIME","OS.TIME.DAYS")][0]
surv=surv[[scol,os_ev,os_t]].dropna()
surv.columns=["sample","OS","OS_time"]; surv=surv.set_index("sample")

# tumor samples only (barcode ...-01)
tum=[s for s in expr.columns if s[-2:] in ("01","11") and s[-2:]=="01"]
E=expr[tum]

def sig_score(genes):
    g=[x for x in genes if x in E.index]
    sub=E.loc[g]                                   # genes x samples
    z=sub.sub(sub.mean(axis=1),axis=0).div(sub.std(axis=1)+1e-9,axis=0)  # z per gene across samples
    return z.mean(axis=0), g                        # Series indexed by sample
latam_score, latam_used = sig_score(LATAM)
df=pd.DataFrame({"LATAM": latam_score})
df["NR4A2"]=E.loc["NR4A2"] if "NR4A2" in E.index else np.nan
df.index=[s[:15] for s in df.index]  # trim to patient-level barcode for join
surv.index=[s[:15] for s in surv.index]
m=df.join(surv, how="inner").dropna(subset=["OS","OS_time","LATAM"])
m=m[m["OS_time"]>0]
print("merged cohort:", m.shape, "| LA-TAM genes used:", latam_used, flush=True)

# correlation LA-TAM vs NR4A2
rho,p=spearmanr(m["LATAM"], m["NR4A2"])
corr=pd.DataFrame([{"spearman_rho":rho,"p":p,"n":len(m)}])
corr.to_csv(f"{TAB}/tcga_kirc_latam_nr4a2_corr.csv", index=False)
print(f"LA-TAM vs NR4A2 spearman rho={rho:.3f} p={p:.2e}", flush=True)

# Cox: continuous LA-TAM and NR4A2 (z), OS
res=[]
for sig in ["LATAM","NR4A2"]:
    d=m[[sig,"OS","OS_time"]].copy(); d[sig]=zscore(d[sig])
    cph=CoxPHFitter().fit(d, duration_col="OS_time", event_col="OS")
    hr=np.exp(cph.params_[sig]); ci=np.exp(cph.confidence_intervals_.loc[sig].values)
    res.append({"signature":sig,"HR_per_SD":hr,"CI_low":ci[0],"CI_high":ci[1],
                "p":cph.summary.loc[sig,"p"],"n":len(d)})
cox=pd.DataFrame(res); cox.to_csv(f"{TAB}/tcga_kirc_survival.csv", index=False)
print(cox.to_string(), flush=True)

# KM by LA-TAM tertile
q=pd.qcut(m["LATAM"],3,labels=["low","mid","high"])
lr=multivariate_logrank_test(m["OS_time"], q, m["OS"])
plt.figure(figsize=(6,5)); kmf=KaplanMeierFitter()
for lev in ["low","mid","high"]:
    idx=q==lev; kmf.fit(m["OS_time"][idx], m["OS"][idx], label=f"LA-TAM {lev} (n={idx.sum()})"); kmf.plot_survival_function()
plt.title(f"TCGA-KIRC OS by LA-TAM signature (logrank p={lr.p_value:.1e})")
plt.xlabel("days"); plt.ylabel("overall survival"); plt.tight_layout()
plt.savefig(f"{FIG}/tcga_kirc_km_latam.png", dpi=150); plt.close()
pd.DataFrame([{"logrank_p_tertile":lr.p_value}]).to_csv(f"{TAB}/tcga_kirc_km_logrank.csv", index=False)
print("TASK6 TCGA DONE", flush=True)
