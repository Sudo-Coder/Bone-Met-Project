#!/usr/bin/env python
"""01_score_clec_lam.py — Phase 1 step 1: define CLEC_LAM CORE state and score it.

Per-object (RCC + prostate myeloid), in each object's own log-normalized space:
  - engineer patient_id / cancer_type / dataset / condition / compartment / is_benign_shared
  - CORE gene detection rates (per cancer, per TAM/TIM)
  - TIM comparability check across cancers (marker profile correlation) BEFORE any TAM+TIM pooling
  - continuous CORE score two ways: sc.tl.score_genes (additive) + decoupler AUCell (rank-based)
  - method concordance (Spearman)
  - thresholded CLEC_LAM-high (top 10/15/20% within TAM+TIM, + 2-comp GMM) — SECONDARY
Writes per-cell score+meta tables for 02_composition_models.py.

Run: envs/rcc_reinterp_venv/bin/python (scanpy + decoupler). Seed 0.
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, decoupler as dc
from scipy.stats import spearmanr
from sklearn.mixture import GaussianMixture

np.random.seed(0); sc.settings.verbosity = 0
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
OUT  = os.path.join(ROOT, "analysis", "rcc_reinterpretation", "outputs")
TAB  = os.path.join(OUT, "tables"); os.makedirs(TAB, exist_ok=True)

CORE = ["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"]
OBJS = {
  "RCC":      dict(path="kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad",
                   label_col="final_label"),
  "prostate": dict(path="prostate-cancer/Cleaned_Data/myeloid_FINAL.h5ad",
                   label_col="scanvi_labels_annotation_model_refined"),
}
# marker panels for the TIM comparability check
MARK = {
  "monocyte":     ["FCN1","VCAN","S100A8","S100A9","LYZ","CD14"],
  "TAM_resident": ["C1QA","C1QB","APOE","MRC1","CD68","MERTK","TREM2"],
  "inflammatory": ["IL1B","CXCL8","CXCL9","CXCL10","TNF","CCL3"],
  "core":         CORE,
}

def parse_patient(sample):
    s = str(sample).split("_",1)[-1].replace(".count","")     # RCC-BM9-Tumor / BMET3-Involved / BMM2-Benign-immune
    parts = s.split("-")
    # patient token = everything up to the site keyword
    site_kw = {"Tumor","Involve","Involved","Noninvolved","Distal","Benign"}
    toks=[]
    for p in parts:
        if p in site_kw or p.startswith("Benign"): break
        toks.append(p)
    pid = "-".join(toks)                                       # RCC-BM9 / BMET3 / BMM2
    pid = pid.replace("RCC-","")                               # BM9
    return pid

def norm_condition(c):
    c=str(c)
    return {"Involve":"Involved","Noninvolved":"Distal"}.get(c, c)

def get_lognorm(a):
    """Return an AnnData on the full gene set with log-normalized X for scoring.
    RCC: raw.X is log-norm (16215 g). prostate: X==raw is log-norm (10809 g)."""
    if a.raw is not None and a.raw.n_vars >= a.n_vars:
        s = a.raw.to_adata()
    else:
        s = a.copy()
    s.obs = a.obs.copy()
    return s

report = {"objects":{}, "core_detection":{}, "tim_comparability":{}, "concordance":{}, "thresholds":{}}
percell = {}

for tag, cfg in OBJS.items():
    a = ad.read_h5ad(os.path.join(ROOT, cfg["path"]))
    lab = cfg["label_col"]
    a.obs["compartment"] = a.obs[lab].astype(str)
    a.obs["condition"]   = a.obs["condition"].astype(str).map(norm_condition)
    a.obs["patient_id"]  = a.obs["Sample"].map(parse_patient)
    a.obs["cancer_type"] = tag
    a.obs["dataset"]     = "rcc_bm" if tag=="RCC" else "prostate_kfoury"
    a.obs["is_benign_shared"] = a.obs["condition"].eq("Benign")
    report["objects"][tag] = dict(n_cells=int(a.n_obs), n_genes=int(a.n_vars),
        compartments=a.obs["compartment"].value_counts().to_dict(),
        conditions=a.obs["condition"].value_counts().to_dict(),
        n_patients=int(a.obs["patient_id"].nunique()),
        patients=sorted(a.obs["patient_id"].unique().tolist()))

    s = get_lognorm(a)
    present = [g for g in CORE if g in s.var_names]
    # ---- CORE detection rates (fraction cells with expr>0) overall + TAM + TIM ----
    det={}
    for grp, mask in [("all", np.ones(s.n_obs,bool)),
                      ("TAM", (a.obs["compartment"]=="TAM").values),
                      ("TIM", (a.obs["compartment"]=="TIM").values)]:
        sub=s[mask]
        X=sub[:,present].X
        X=X.toarray() if hasattr(X,"toarray") else np.asarray(X)
        det[grp]={g:float((X[:,i]>0).mean()) for i,g in enumerate(present)}
        det[grp]["_n"]=int(mask.sum())
    report["core_detection"][tag]=det

    # ---- TIM comparability marker profile (mean log-norm expr per compartment) ----
    prof={}
    for panel,genes in MARK.items():
        gg=[g for g in genes if g in s.var_names]
        for comp in ["TAM","TIM"]:
            m=(a.obs["compartment"]==comp).values
            if m.sum()==0: continue
            X=s[m][:,gg].X; X=X.toarray() if hasattr(X,"toarray") else np.asarray(X)
            prof[f"{comp}|{panel}"]=float(X.mean())
    report["tim_comparability"][tag]={"profile_means":prof}
    # store TIM per-gene mean vector for cross-cancer correlation
    tim_genes=[g for g in (MARK["monocyte"]+MARK["TAM_resident"]+MARK["inflammatory"]) if g in s.var_names]
    m=(a.obs["compartment"]=="TIM").values
    X=s[m][:,tim_genes].X; X=X.toarray() if hasattr(X,"toarray") else np.asarray(X)
    report["tim_comparability"][tag]["tim_gene_means"]=dict(zip(tim_genes, X.mean(0).tolist()))

    # ---- scoring ----
    sc.tl.score_genes(s, present, score_name="CLEC_LAM_addmodule", random_state=0)
    net=pd.DataFrame({"source":"CLEC_LAM_CORE","target":present,"weight":1.0})
    dc.mt.aucell(s, net, tmin=3, verbose=False)
    s.obs["CLEC_LAM_aucell"]=s.obsm["score_aucell"]["CLEC_LAM_CORE"].values

    rho,_=spearmanr(s.obs["CLEC_LAM_addmodule"], s.obs["CLEC_LAM_aucell"])
    # concordance within TAM+TIM too
    tt=a.obs["compartment"].isin(["TAM","TIM"]).values
    rho_tt,_=spearmanr(s.obs["CLEC_LAM_addmodule"][tt], s.obs["CLEC_LAM_aucell"][tt])
    report["concordance"][tag]=dict(spearman_all=float(rho), spearman_TAM_TIM=float(rho_tt), n_TAM_TIM=int(tt.sum()))

    # ---- thresholded high-cells within TAM+TIM (SECONDARY) ----
    thr={}
    score=s.obs["CLEC_LAM_aucell"].values
    idx_tt=np.where(tt)[0]
    for q in [0.90,0.85,0.80]:
        pct=int(round((1-q)*100))
        cut=np.quantile(score[idx_tt], q)
        s.obs[f"high_top{pct}"]=(score>=cut)&tt
        thr[f"top{pct}pct"]=dict(cut=float(cut), n_high=int(((score>=cut)&tt).sum()))
    gm=GaussianMixture(2,random_state=0).fit(score[idx_tt].reshape(-1,1))
    hi_comp=int(np.argmax(gm.means_.ravel()))
    lab_gm=gm.predict(score.reshape(-1,1))
    s.obs["high_gmm"]=(lab_gm==hi_comp)&tt
    thr["gmm"]=dict(n_high=int(s.obs["high_gmm"].sum()),
                    means=gm.means_.ravel().tolist(), weights=gm.weights_.tolist())
    report["thresholds"][tag]=thr

    # ---- save per-cell table ----
    keep=["patient_id","cancer_type","dataset","condition","compartment","is_benign_shared","Sample",
          "CLEC_LAM_addmodule","CLEC_LAM_aucell","high_top10","high_top15","high_top20","high_gmm"]
    df=s.obs[keep].copy(); df.index.name="cell"
    percell[tag]=df

pc=pd.concat(percell.values())
pc.to_parquet(os.path.join(TAB,"clec_lam_percell.parquet"))

# ---- cross-cancer TIM comparability correlation ----
r_g=set(report["tim_comparability"]["RCC"]["tim_gene_means"])
p_g=set(report["tim_comparability"]["prostate"]["tim_gene_means"])
shared=sorted(r_g & p_g)
rv=[report["tim_comparability"]["RCC"]["tim_gene_means"][g] for g in shared]
pv=[report["tim_comparability"]["prostate"]["tim_gene_means"][g] for g in shared]
tim_rho,_=spearmanr(rv,pv)
report["tim_comparability"]["cross_cancer_TIM_marker_spearman"]=dict(rho=float(tim_rho), n_genes=len(shared), genes=shared)

with open(os.path.join(TAB,"phase1_scoring_report.json"),"w") as f:
    json.dump(report,f,indent=2)

print("=== CORE detection (fraction>0) ===")
for tag in OBJS:
    d=report["core_detection"][tag]
    print(f"[{tag}] TAM n={d['TAM']['_n']} TIM n={d['TIM']['_n']}")
    for grp in ["all","TAM","TIM"]:
        print(f"  {grp}:", {g:round(d[grp][g],2) for g in CORE if g in d[grp]})
print("\n=== method concordance (Spearman addmodule vs AUCell) ===")
for tag in OBJS: print(f"  {tag}: all={report['concordance'][tag]['spearman_all']:.3f}  TAM+TIM={report['concordance'][tag]['spearman_TAM_TIM']:.3f}")
print(f"\n=== cross-cancer TIM marker-profile Spearman: rho={tim_rho:.3f} (n={len(shared)} genes) ===")
print("\nwrote:", os.path.join(TAB,"clec_lam_percell.parquet"), "+ phase1_scoring_report.json")
