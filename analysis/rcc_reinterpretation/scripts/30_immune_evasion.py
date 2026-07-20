#!/usr/bin/env python
import os, warnings, json
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, decoupler as dc
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr
np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables")
A=ad.read_h5ad(os.path.join(ROOT,"kidney-cancer/Cleaned_Data/integrated.h5ad"))
s=A.raw.to_adata() if (A.raw is not None and A.raw.n_vars>=A.n_vars) else A.copy(); s.obs=A.obs.copy()

MOD={"complement_C1Q":["C1QA","C1QB","C1QC"],"complement_C1Q_C3":["C1QA","C1QB","C1QC","C3"],
     "CLEC_LAM8":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"],
     "RCC_skew_CORE":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2"],
     "APOE_TREM2":["APOE","APOC1","TREM2","TYROBP"],"MERTK_GPNMB":["MERTK","GPNMB"],
     "SPP1_TAM":["SPP1","FN1","MMP9","CTSB","CTSD"],"inflammatory_mono":["IL1B","CXCL8","S100A8","S100A9","FCN1","VCAN"],
     "panTAM":["CD68","CD163","MRC1","CSF1R","LYZ"]}
OUT={"CD8_exhaustion":["TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1"],
     "cytotoxicity":["GZMB","PRF1","NKG7","GNLY","IFNG"],
     "MHC_I_APM":["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","NLRC5"],
     "MHC_II_APC":["HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","CIITA"]}
net=pd.concat([pd.DataFrame({"source":k,"target":[g for g in v if g in s.var_names],"weight":1.0}) for k,v in {**MOD,**OUT}.items()])
dc.mt.aucell(s,net,tmin=2,verbose=False); SC=s.obsm["score_aucell"]
ct=A.obs["cell type"].astype(str).values
cond={"Involve":"Involved","Noninvolved":"Distal"}
A.obs["condition"]=A.obs["condition"].astype(str).map(lambda c:cond.get(c,c))
def pid(x):
    x=str(x).split("_",1)[-1].replace(".count","")
    for k in ["-Tumor","-Involve","-Noninvolved","-Distal","-Benign"]: x=x.split(k)[0]
    return x.replace("RCC-","")
A.obs["pid"]=A.obs["Sample"].map(pid)
MYELOID={"Mono1","Mono2","Mono3","TAM","TIM","mDC","Monocyte Pro","Mono-1","Mono-2","Mono-3"}
CD8={"CTL-1","CTL-2","CTL-3","CD8 Naive"}; TCELL={"Thelper","CD4 Naive","CD8 Naive","CTL-1","CTL-2","CTL-3","NKT","Treg","Proliferating T"}
APC=MYELOID|{"Mature B","Immature B cells","memBcell","Pro Bcell","PDC","mDC"}

df=pd.DataFrame(index=A.obs_names);
for c in list(MOD)+list(OUT): df[c]=SC[c].values if c in SC.columns else np.nan
df["ct"]=ct; df["pid"]=A.obs["pid"].values; df["cond"]=A.obs["condition"].values; df["Sample"]=A.obs["Sample"].values
def agg(sub):
    d={}
    my=sub[sub.ct.isin(MYELOID)];
    for m in MOD: d[m]=my[m].mean() if len(my) else np.nan
    cd8=sub[sub.ct.isin(CD8)]; d["CD8_exhaustion"]=cd8["CD8_exhaustion"].mean() if len(cd8)>=5 else np.nan
    d["cytotoxicity"]=cd8["cytotoxicity"].mean() if len(cd8)>=5 else np.nan
    tum=sub[sub.ct=="Tumor"]; d["MHC_I_APM"]=tum["MHC_I_APM"].mean() if len(tum)>=5 else np.nan
    apc=sub[sub.ct.isin(APC)]; d["MHC_II_APC"]=apc["MHC_II_APC"].mean() if len(apc)>=5 else np.nan
    n=len(sub); t=sub[sub.ct.isin(TCELL)]
    d["Treg_fraction"]=(sub.ct=="Treg").sum()/max(1,len(t))
    d["TAM_fraction"]=sub.ct.isin(["TAM","TIM"]).sum()/n
    d["CD8_fraction"]=sub.ct.isin(CD8).sum()/n
    d["malignant_fraction"]=(sub.ct=="Tumor").sum()/n
    d["n_cells"]=n
    return pd.Series(d)
S=df.groupby(["Sample","pid","cond"],observed=True).apply(agg).reset_index()
S.to_csv(os.path.join(TAB,"phase3_sample_table.csv"),index=False)

CTRLS=["panTAM","SPP1_TAM","inflammatory_mono","MERTK_GPNMB","TAM_fraction"]
TESTMOD=["complement_C1Q","complement_C1Q_C3","CLEC_LAM8","RCC_skew_CORE","APOE_TREM2","MERTK_GPNMB"]+["panTAM","SPP1_TAM","inflammatory_mono"]
OUTCOMES=["CD8_exhaustion","cytotoxicity","Treg_fraction","MHC_I_APM","MHC_II_APC"]
def zscore(x):
    x=pd.to_numeric(x,errors="coerce"); return (x-x.mean())/x.std()
rows=[]
for oc in OUTCOMES:
    D=S.dropna(subset=[oc]).copy()
    if len(D)<8:
        for m in TESTMOD: rows.append(dict(outcome=oc,module=m,coef=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(D),note="too few"))
        continue
    for cvar in ["TAM_fraction","CD8_fraction","malignant_fraction",oc]+TESTMOD:
        if cvar in D: D[cvar+"_z"] if False else None
    for m in TESTMOD:
        d=D.dropna(subset=[m,"TAM_fraction","CD8_fraction","malignant_fraction"]).copy()
        if len(d)<8: rows.append(dict(outcome=oc,module=m,coef=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d),note="too few")); continue
        for c in [m,"TAM_fraction","CD8_fraction","malignant_fraction"]: d[c+"_z"]=zscore(d[c])
        d["y"]=zscore(d[oc])
        try:
            fit=smf.ols(f"y ~ {m}_z + TAM_fraction_z + CD8_fraction_z + malignant_fraction_z + C(cond)",d).fit()
            ci=fit.conf_int().loc[f"{m}_z"]
            rows.append(dict(outcome=oc,module=m,coef=fit.params[f"{m}_z"],ci_low=ci[0],ci_high=ci[1],p=fit.pvalues[f"{m}_z"],n=len(d),note="adj z"))
        except Exception as e:
            rows.append(dict(outcome=oc,module=m,coef=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d),note=str(e)[:30]))
R=pd.DataFrame(rows)

R["fdr"]=np.nan
for oc in OUTCOMES:
    mask=(R.outcome==oc)&R.p.notna()
    if mask.sum(): R.loc[mask,"fdr"]=multipletests(R.loc[mask,"p"],method="fdr_bh")[1]
R.to_csv(os.path.join(TAB,"phase3_adjusted_associations.csv"),index=False)

h2h=[]
for oc in OUTCOMES:
    d=S.dropna(subset=[oc,"complement_C1Q","panTAM","TAM_fraction","CD8_fraction","malignant_fraction"]).copy()
    if len(d)<8: continue
    for c in ["complement_C1Q","panTAM","TAM_fraction","CD8_fraction","malignant_fraction"]: d[c+"_z"]=zscore(d[c])
    d["y"]=zscore(d[oc])
    try:
        fit=smf.ols("y ~ complement_C1Q_z + panTAM_z + TAM_fraction_z + CD8_fraction_z + malignant_fraction_z + C(cond)",d).fit()
        h2h.append(dict(outcome=oc,n=len(d),
            complement_coef=fit.params["complement_C1Q_z"],complement_p=fit.pvalues["complement_C1Q_z"],
            panTAM_coef=fit.params["panTAM_z"],panTAM_p=fit.pvalues["panTAM_z"]))
    except Exception: pass
H=pd.DataFrame(h2h); H.to_csv(os.path.join(TAB,"phase3_headtohead_complement_vs_panTAM.csv"),index=False)

pd.set_option("display.width",220,"display.max_colwidth",18)
print("=== sample table (RCC) n samples:",len(S),"| by condition:",dict(S.cond.value_counts()),"===")
print("\n=== adjusted associations (module coef, z; outcome ~ module + TAMfrac + CD8frac + malignfrac + cond) ===")
print(R[["outcome","module","coef","ci_low","ci_high","p","fdr","n"]].round(3).to_string(index=False))
print("\n=== head-to-head complement vs panTAM (same model) ===")
print(H.round(3).to_string(index=False))
