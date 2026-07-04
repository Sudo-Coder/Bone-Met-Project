#!/usr/bin/env python
"""31_pca_immune_evasion.py — PCa mirror of Phase 3 immune-evasion (prostate full-niche).
Per prostate sample: myeloid module scores vs CD8 exhaustion/cytotoxicity, Treg fraction, MHC-I/APM,
MHC-II/APC; adjusted for TAM/CD8/malignant fraction + condition. Contrast-primary (complement vs panTAM
head-to-head), mirroring RCC 30_immune_evasion. Run: envs/rcc_reinterp_venv/bin/python. Seed 0.
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, decoupler as dc
import statsmodels.formula.api as smf
np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","pca_comparator","outputs","tables"); os.makedirs(TAB,exist_ok=True)
A=ad.read_h5ad(os.path.join(ROOT,"prostate-cancer/Cleaned_Data/integrated_with_kfoury_labels.h5ad"))
s=A.raw.to_adata() if (A.raw is not None and A.raw.n_vars>=A.n_vars) else A.copy(); s.obs=A.obs.copy()
MOD={"complement_C1Q":["C1QA","C1QB","C1QC"],"CLEC_LAM8":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"],
     "RCC_skew_CORE":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2"],"APOE_TREM2":["APOE","APOC1","TREM2","TYROBP"],
     "MERTK_GPNMB":["MERTK","GPNMB"],"SPP1_TAM":["SPP1","FN1","MMP9","CTSB","CTSD"],"panTAM":["CD68","CD163","MRC1","CSF1R","LYZ"],
     "ATF3_NFkB":["ATF3","NFKB1","RELA","JUN","JUNB","DUSP1","TNFAIP3","NFKBIA"]}
OUT={"CD8_exhaustion":["TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1"],"cytotoxicity":["GZMB","PRF1","NKG7","GNLY","IFNG"],
     "MHC_I_APM":["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","NLRC5"],"MHC_II_APC":["HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","CIITA"]}
net=pd.concat([pd.DataFrame({"source":k,"target":[g for g in v if g in s.var_names],"weight":1.0}) for k,v in {**MOD,**OUT}.items()])
dc.mt.aucell(s,net,tmin=2,verbose=False); SC=s.obsm["score_aucell"]
ct=A.obs["kfoury_annotation"].astype(str).values
def cond(x):
    x=str(x)
    if "Benign" in x: return "Benign"
    if "-Tumor" in x: return "Tumor"
    if "-Involved" in x: return "Involved"
    if "-Distal" in x: return "Distal"
    return "other"
def pid(x):
    x=str(x).split("_",1)[-1].replace(".count","")
    for k in ["-Tumor","-Involved","-Distal","-Benign"]: x=x.split(k)[0]
    return x
MYE={"Mono1","Mono2","Mono3","TAM","TIM","mDC","Monocyte prog","pDC"}
CD8={"CTL-1","CTL-2","CD8+ Naive"}; TREG={"Treg Active","Treg Resting"}
TCELL={"CD4+ Naive","CTL-1","CTL-2","CD8+ Naive","Th1/17","NKT","NK","Treg Active","Treg Resting"}
APC=MYE|{"Mature B","Immature B cells","memBcell","Pro-B"}
df=pd.DataFrame(index=A.obs_names)
for c in list(MOD)+list(OUT): df[c]=SC[c].values if c in SC.columns else np.nan
df["ct"]=ct; df["cond"]=A.obs["Sample"].map(cond).values; df["pid"]=A.obs["Sample"].map(pid).values; df["Sample"]=A.obs["Sample"].values
def agg(sub):
    d={}; my=sub[sub.ct.isin(MYE)]
    for m in MOD: d[m]=my[m].mean() if len(my) else np.nan
    cd8=sub[sub.ct.isin(CD8)]; d["CD8_exhaustion"]=cd8["CD8_exhaustion"].mean() if len(cd8)>=5 else np.nan
    d["cytotoxicity"]=cd8["cytotoxicity"].mean() if len(cd8)>=5 else np.nan
    tum=sub[sub.ct=="Tumor"]; d["MHC_I_APM"]=tum["MHC_I_APM"].mean() if len(tum)>=5 else np.nan
    apc=sub[sub.ct.isin(APC)]; d["MHC_II_APC"]=apc["MHC_II_APC"].mean() if len(apc)>=5 else np.nan
    n=len(sub); t=sub[sub.ct.isin(TCELL)]
    d["Treg_fraction"]=sub.ct.isin(TREG).sum()/max(1,len(t))
    d["TAM_fraction"]=sub.ct.isin(["TAM","TIM"]).sum()/n; d["CD8_fraction"]=sub.ct.isin(CD8).sum()/n; d["malignant_fraction"]=(sub.ct=="Tumor").sum()/n; d["n_cells"]=n
    return pd.Series(d)
S=df.groupby(["Sample","pid","cond"],observed=True).apply(agg).reset_index(); S.to_csv(os.path.join(TAB,"pca_phase3_sample_table.csv"),index=False)
def zc(x): x=pd.to_numeric(x,errors="coerce"); return (x-x.mean())/x.std()
TEST=["complement_C1Q","CLEC_LAM8","RCC_skew_CORE","APOE_TREM2","MERTK_GPNMB","panTAM","SPP1_TAM","ATF3_NFkB"]
OUTC=["CD8_exhaustion","cytotoxicity","Treg_fraction","MHC_I_APM","MHC_II_APC"]
rows=[]
for oc in OUTC:
    D=S.dropna(subset=[oc]).copy()
    for m in TEST:
        d=D.dropna(subset=[m,"TAM_fraction","CD8_fraction","malignant_fraction"]).copy()
        if len(d)<8: rows.append(dict(outcome=oc,module=m,coef=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d))); continue
        for c in [m,"TAM_fraction","CD8_fraction","malignant_fraction"]: d[c+"_z"]=zc(d[c])
        d["y"]=zc(d[oc])
        try:
            fit=smf.ols(f"y ~ {m}_z + TAM_fraction_z + CD8_fraction_z + malignant_fraction_z + C(cond)",d).fit()
            ci=fit.conf_int().loc[f"{m}_z"]; rows.append(dict(outcome=oc,module=m,coef=fit.params[f"{m}_z"],ci_low=ci[0],ci_high=ci[1],p=fit.pvalues[f"{m}_z"],n=len(d)))
        except Exception as e: rows.append(dict(outcome=oc,module=m,coef=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d)))
R=pd.DataFrame(rows); R.to_csv(os.path.join(TAB,"pca_phase3_adjusted_associations.csv"),index=False)
h2h=[]
for oc in OUTC:
    d=S.dropna(subset=[oc,"complement_C1Q","panTAM","TAM_fraction","CD8_fraction","malignant_fraction"]).copy()
    if len(d)<8: continue
    for c in ["complement_C1Q","panTAM","TAM_fraction","CD8_fraction","malignant_fraction"]: d[c+"_z"]=zc(d[c])
    d["y"]=zc(d[oc])
    try:
        fit=smf.ols("y ~ complement_C1Q_z + panTAM_z + TAM_fraction_z + CD8_fraction_z + malignant_fraction_z + C(cond)",d).fit()
        h2h.append(dict(outcome=oc,n=len(d),complement_coef=fit.params["complement_C1Q_z"],complement_p=fit.pvalues["complement_C1Q_z"],panTAM_coef=fit.params["panTAM_z"],panTAM_p=fit.pvalues["panTAM_z"]))
    except Exception: pass
H=pd.DataFrame(h2h); H.to_csv(os.path.join(TAB,"pca_phase3_headtohead.csv"),index=False)
pd.set_option("display.width",220,"display.max_colwidth",16)
print("PCa full-niche samples:",len(S),"| by condition:",dict(S.cond.value_counts()))
print("\n=== PCa adjusted immune-evasion associations ===")
print(R[["outcome","module","coef","ci_low","ci_high","p","n"]].round(3).to_string(index=False))
print("\n=== PCa head-to-head complement vs panTAM ===")
print(H.round(3).to_string(index=False))
