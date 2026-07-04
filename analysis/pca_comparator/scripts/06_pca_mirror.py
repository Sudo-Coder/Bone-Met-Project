#!/usr/bin/env python
"""06_pca_mirror.py — LOCKED comparator: is the RCC complement_C1Q program present in PCa myeloid?
Mirrors the RCC framework (04/02) for prostate. Scores modules (+ATF3/NFkB) on both myeloid objects,
patient-level pseudobulk; runs PCa tumor-vs-benign + tumor/involved/distal gradient + within-PCa-tumor
TAM/TIM localization of complement_C1Q; emits a PCa-vs-RCC module contrast table.
Run: envs/rcc_reinterp_venv/bin/python. Seed 0.
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, decoupler as dc
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","pca_comparator","outputs","tables"); os.makedirs(TAB,exist_ok=True); MIN=20
MODULES={
 "CLEC_LAM8":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"],
 "complement_C1Q":["C1QA","C1QB","C1QC"],
 "RCC_skew_CORE":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2"],
 "APOE_TREM2":["APOE","APOC1","TREM2","TYROBP"],
 "MERTK_GPNMB":["MERTK","GPNMB"],
 "SPP1_TAM":["SPP1","FN1","MMP9","CTSB","CTSD"],
 "panTAM":["CD68","CD163","MRC1","CSF1R","LYZ"],
 "ATF3_IEG":["ATF3","JUN","JUNB","JUND","FOS","FOSB","EGR1","NR4A1","DUSP1"],
 "NFkB":["NFKB1","RELA","RELB","NFKBIA","NFKBIZ","TNFAIP3","CCL3","CCL4"],
 "ATF3_NFkB":["ATF3","NFKB1","RELA","JUN","JUNB","DUSP1","TNFAIP3","NFKBIA"],
}
OBJS={"RCC":dict(path="kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad",lab="final_label"),
      "prostate":dict(path="prostate-cancer/Cleaned_Data/myeloid_FINAL.h5ad",lab="scanvi_labels_annotation_model_refined")}
def ppid(s):
    s=str(s).split("_",1)[-1].replace(".count","")
    for k in ["-Tumor","-Involve","-Noninvolved","-Distal","-Benign"]: s=s.split(k)[0]
    return s.replace("RCC-","")
def nc(c): return {"Involve":"Involved","Noninvolved":"Distal"}.get(str(c),str(c))
per=[]
for tag,cfg in OBJS.items():
    a=ad.read_h5ad(os.path.join(ROOT,cfg["path"]))
    s=a.raw.to_adata() if (a.raw is not None and a.raw.n_vars>=a.n_vars) else a.copy(); s.obs=a.obs.copy()
    net=pd.concat([pd.DataFrame({"source":m,"target":[g for g in v if g in s.var_names],"weight":1.0}) for m,v in MODULES.items()])
    dc.mt.aucell(s,net,tmin=2,verbose=False); SC=s.obsm["score_aucell"]
    df=pd.DataFrame(index=a.obs_names)
    for m in MODULES: df[m]=SC[m].values if m in SC.columns else np.nan
    df["compartment"]=a.obs[cfg["lab"]].astype(str).values
    df["condition"]=a.obs["condition"].astype(str).map(nc).values
    df["patient_id"]=a.obs["Sample"].map(ppid).values; df["cancer_type"]=tag; df["Sample"]=a.obs["Sample"].values
    per.append(df)
pc=pd.concat(per); pc.to_parquet(os.path.join(TAB,"pca_mirror_module_scores_percell.parquet"))

def pseudobulk(df):
    g=df.groupby(["Sample","cancer_type","patient_id","condition"],observed=True)
    return g.agg(**{**{m:(m,"mean") for m in MODULES},"n_cells":("condition","size")}).reset_index()
pb=pseudobulk(pc); pb=pb[pb.n_cells>=MIN]
pb_tt=pseudobulk(pc[pc.compartment.isin(["TAM","TIM"])]); pb_tt=pb_tt[pb_tt.n_cells>=10]
pb.to_csv(os.path.join(TAB,"pca_mirror_pseudobulk_all_myeloid.csv"),index=False)

def tvb(dfc,m):
    """tumor-vs-benign Delta within one cancer (all-myeloid pseudobulk), mixed +patient RE; SD units."""
    d=dfc[dfc.condition.isin(["Tumor","Benign"])].dropna(subset=[m]).copy()
    if d[m].nunique()<3 or (d.condition=="Benign").sum()<3 or (d.condition=="Tumor").sum()<3: return (np.nan,)*5
    d["condition"]=pd.Categorical(d["condition"],["Benign","Tumor"]); sd=float(d[m].std())
    try:
        mm=smf.mixedlm(f"{m} ~ C(condition, Treatment('Benign'))",d,groups=d["patient_id"]).fit(reml=False,method="lbfgs")
        e,se,p=mm.params["C(condition, Treatment('Benign'))[T.Tumor]"],mm.bse["C(condition, Treatment('Benign'))[T.Tumor]"],mm.pvalues["C(condition, Treatment('Benign'))[T.Tumor]"]
    except Exception:
        mm=smf.ols(f"{m} ~ C(condition, Treatment('Benign'))",d).fit(cov_type="cluster",cov_kwds={"groups":d["patient_id"]})
        ci=mm.conf_int().loc["C(condition, Treatment('Benign'))[T.Tumor]"]; e,se,p=mm.params[1],(ci[1]-ci[0])/3.92,mm.pvalues[1]
    return e,e-1.96*se,e+1.96*se,p,sd
def gradient(dfc,m):
    """linear trend across ordered condition Benign<Distal<Involved<Tumor (all-myeloid), patient cluster."""
    order={"Benign":0,"Distal":1,"Involved":2,"Tumor":3}
    d=dfc[dfc.condition.isin(order)].dropna(subset=[m]).copy(); d["ord"]=d.condition.map(order)
    if len(d)<10: return (np.nan,np.nan,np.nan)
    fit=smf.ols(f"{m} ~ ord",d).fit(cov_type="cluster",cov_kwds={"groups":d["patient_id"]})
    ci=fit.conf_int().loc["ord"]; return fit.params["ord"],fit.pvalues["ord"],(ci[0],ci[1])

rows=[]
INT_form="{m} ~ C(cancer_type, Treatment('prostate'))*C(condition, Treatment('Benign'))"
INT_coef="C(cancer_type, Treatment('prostate'))[T.RCC]:C(condition, Treatment('Benign'))[T.Tumor]"
for m in MODULES:
    rcc=tvb(pb[pb.cancer_type=="RCC"],m); pca=tvb(pb[pb.cancer_type=="prostate"],m)
    gr=gradient(pb[pb.cancer_type=="prostate"],m)
    # interaction
    dd=pb[pb.condition.isin(["Tumor","Benign"])].dropna(subset=[m]).copy()
    dd["condition"]=pd.Categorical(dd["condition"],["Benign","Tumor"]); dd["cancer_type"]=pd.Categorical(dd["cancer_type"],["prostate","RCC"])
    try:
        mm=smf.mixedlm(INT_form.format(m=m),dd,groups=dd["patient_id"]).fit(reml=False,method="lbfgs"); ie,ip=mm.params[INT_coef],mm.pvalues[INT_coef]
    except Exception:
        mm=smf.ols(INT_form.format(m=m),dd).fit(cov_type="cluster",cov_kwds={"groups":dd["patient_id"]}); ie,ip=mm.params[INT_coef],mm.pvalues[INT_coef]
    rows.append(dict(module=m,
        RCC_TvB=rcc[0],RCC_TvB_SD=rcc[0]/rcc[4] if rcc[4] else np.nan,RCC_p=rcc[3],
        PCa_TvB=pca[0],PCa_TvB_ci=f"[{pca[1]:.3f},{pca[2]:.3f}]" if pca[0]==pca[0] else "",PCa_TvB_SD=pca[0]/pca[4] if pca[4] else np.nan,PCa_p=pca[3],
        PCa_gradient_slope=gr[0],PCa_gradient_p=gr[1],
        interaction_RCCvsPCa=ie,interaction_p=ip))
C=pd.DataFrame(rows)
C["PCa_TvB_FDR"]=multipletests(C["PCa_p"].fillna(1),method="fdr_bh")[1]
C.to_csv(os.path.join(TAB,"pca_vs_rcc_module_contrast.csv"),index=False)

# within-PCa-tumor localization: complement_C1Q (+key) in TAM+TIM vs Mono (tumor samples)
loc=[]
pt=pc[(pc.cancer_type=="prostate")&(pc.condition=="Tumor")].copy()
pt["grp"]=np.where(pt.compartment.isin(["TAM","TIM"]),"TAM_TIM",np.where(pt.compartment.isin(["Mono1","Mono2","Mono3"]),"Mono","other"))
for m in ["complement_C1Q","CLEC_LAM8","APOE_TREM2","SPP1_TAM","ATF3_NFkB"]:
    # patient-level mean per grp
    g=pt[pt.grp.isin(["TAM_TIM","Mono"])].groupby(["patient_id","grp"],observed=True)[m].mean().reset_index()
    piv=g.pivot(index="patient_id",columns="grp",values=m).dropna()
    if len(piv)>=4:
        from scipy.stats import wilcoxon
        try: _,p=wilcoxon(piv["TAM_TIM"],piv["Mono"])
        except Exception: p=np.nan
        loc.append(dict(module=m,n_patients=len(piv),TAM_TIM_mean=float(piv["TAM_TIM"].mean()),Mono_mean=float(piv["Mono"].mean()),delta=float((piv["TAM_TIM"]-piv["Mono"]).mean()),wilcoxon_p=p))
L=pd.DataFrame(loc); L.to_csv(os.path.join(TAB,"pca_within_tumor_localization.csv"),index=False)

pd.set_option("display.width",240,"display.max_colwidth",20)
print("=== PCa tumor-vs-benign + PCa-vs-RCC contrast (all-myeloid pseudobulk) ===")
print(C[["module","RCC_TvB_SD","RCC_p","PCa_TvB","PCa_TvB_SD","PCa_p","PCa_TvB_FDR","PCa_gradient_slope","PCa_gradient_p","interaction_RCCvsPCa","interaction_p"]].round(3).to_string(index=False))
print("\n=== within-PCa-tumor localization: TAM+TIM vs Mono (patient-level) ===")
print(L.round(3).to_string(index=False))
