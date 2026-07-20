#!/usr/bin/env python
import os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, decoupler as dc
import statsmodels.formula.api as smf, statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables"); os.makedirs(TAB,exist_ok=True)
CORE=["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"]
MIN_CELLS=20

pc=pd.read_parquet(os.path.join(TAB,"clec_lam_percell.parquet"))

for ct,sub in pc.groupby("cancer_type",observed=True):
    for pct in (10,15,20):
        cut=np.quantile(sub["CLEC_LAM_aucell"], 1-pct/100)
        pc.loc[sub.index, f"hi{pct}"]=(sub["CLEC_LAM_aucell"]>=cut)
for pct in (10,15,20): pc[f"hi{pct}"]=pc[f"hi{pct}"].astype(bool)

tt=pc[pc.compartment.isin(["TAM","TIM"])]
avail=(tt.groupby(["cancer_type","condition","Sample"],observed=True).size()
         .groupby(level=[0,1]).agg(n_samples="size", n_ge10=lambda x:(x>=10).sum(),
                                   median_cells="median")).reset_index()
avail.to_csv(os.path.join(TAB,"tamtim_availability.csv"),index=False)

def pseudobulk(df):
    g=df.groupby(["Sample","cancer_type","dataset","patient_id","condition"],observed=True)
    return g.agg(score=("CLEC_LAM_aucell","mean"), score_add=("CLEC_LAM_addmodule","mean"),
                 n_cells=("CLEC_LAM_aucell","size"),
                 n_hi10=("hi10","sum"), n_hi15=("hi15","sum"), n_hi20=("hi20","sum")).reset_index()
pb=pseudobulk(pc); pb["unit"]="all_myeloid"; pb["ok"]=pb["n_cells"]>=MIN_CELLS
pb_tt=pseudobulk(tt); pb_tt["unit"]="TAM+TIM"; pb_tt["ok"]=pb_tt["n_cells"]>=10
pb.to_csv(os.path.join(TAB,"pseudobulk_all_myeloid.csv"),index=False)
pb_tt.to_csv(os.path.join(TAB,"pseudobulk_TAMTIM.csv"),index=False)

results=[]
def add(name,tier,est,lo,hi,p,extra=""):
    results.append(dict(test=name,tier=tier,estimate=est,ci_low=lo,ci_high=hi,p=p,note=extra))
def mixedlm_coef(df,formula,group,coef):
    m=smf.mixedlm(formula,df,groups=df[group]).fit(reml=False,method="lbfgs")
    return m.params[coef], m.params[coef]-1.96*m.bse[coef], m.params[coef]+1.96*m.bse[coef], m.pvalues[coef]
def ols_coef(df,formula,coef,cluster=None):
    m=(smf.ols(formula,df).fit(cov_type="cluster",cov_kwds={"groups":df[cluster]}) if cluster
       else smf.ols(formula,df).fit())
    ci=m.conf_int().loc[coef]; return m.params[coef],ci[0],ci[1],m.pvalues[coef]

form="score ~ C(cancer_type, Treatment('prostate'))*C(condition, Treatment('Benign'))"
coef="C(cancer_type, Treatment('prostate'))[T.RCC]:C(condition, Treatment('Benign'))[T.Tumor]"

rcc=pb[(pb.cancer_type=="RCC")&pb.ok].copy()
rcc["condition"]=pd.Categorical(rcc["condition"],["Benign","Distal","Involved","Tumor"])
cf="score ~ C(condition, Treatment('Benign'))"; tgt="C(condition, Treatment('Benign'))[T.Tumor]"
try: e,l,h,p=mixedlm_coef(rcc,cf,"patient_id",tgt)
except Exception: e,l,h,p=ols_coef(rcc,cf,tgt,cluster="patient_id")
tm=rcc.loc[rcc.condition=="Tumor","score"]; bn=rcc.loc[rcc.condition=="Benign","score"]
add("C1_RCC_tumor_vs_benign_allmyeloid","tumor-enriched (patient-level)",e,l,h,p,
    f"n_tumor={len(tm)} n_benign={len(bn)} tumor_mean={tm.mean():.4f} benign_mean={bn.mean():.4f}")

b=pb[pb.condition.isin(["Tumor","Benign"])&pb.ok].copy()
b["condition"]=pd.Categorical(b["condition"],["Benign","Tumor"])
b["cancer_type"]=pd.Categorical(b["cancer_type"],["prostate","RCC"])
try: e,l,h,p=mixedlm_coef(b,form,"patient_id",coef); meth="mixedlm"
except Exception: e,l,h,p=ols_coef(b,form,coef,cluster="patient_id"); meth="ols_cluster"
add("C2_interaction_ownBenign_allmyeloid","RCC-skewed (interaction, HEADLINE)",e,l,h,p,f"method={meth}")
e2,l2,h2,p2=ols_coef(b,form,coef,cluster="patient_id")
add("C2b_interaction_ownBenign_OLScluster","RCC-skewed (sensitivity)",e2,l2,h2,p2)

loo=[]
for pid in b.patient_id.unique():
    try:
        ee,_,_,pp=ols_coef(b[b.patient_id!=pid],form,coef,cluster="patient_id"); loo.append((pid,ee,pp))
    except Exception: pass
loo=pd.DataFrame(loo,columns=["dropped_patient","interaction_est","p"]); loo.to_csv(os.path.join(TAB,"C2_leave_one_out.csv"),index=False)

fr=b.copy(); fr["succ"]=fr["n_hi15"].astype(int); fr["fail"]=(fr["n_cells"]-fr["n_hi15"]).astype(int)
fr["frac"]=fr["n_hi15"]/fr["n_cells"]; fr["asin"]=np.arcsin(np.sqrt(fr["frac"].clip(0,1)))
try:
    gb=smf.glm("succ + fail ~ C(cancer_type, Treatment('prostate'))*C(condition, Treatment('Benign'))",
               data=fr,family=sm.families.Binomial()).fit(cov_type="cluster",cov_kwds={"groups":fr["patient_id"]})
    ci=gb.conf_int().loc[coef]
    add("C3_fraction_interaction_binomialGLM","RCC-skewed fraction (interaction)",
        float(gb.params[coef]),float(ci[0]),float(ci[1]),float(gb.pvalues[coef]),"log-OR CLEC_LAM-hi; cluster SE")
except Exception as ex:
    add("C3_fraction_interaction_binomialGLM","RCC-skewed fraction",np.nan,np.nan,np.nan,np.nan,f"failed:{ex}")
ea,la,ha,pa=ols_coef(fr,"asin ~ C(cancer_type, Treatment('prostate'))*C(condition, Treatment('Benign'))",coef,cluster="patient_id")
add("C3b_fraction_interaction_arcsinLMM","RCC-skewed fraction (arcsin sensitivity)",ea,la,ha,pa)

import scipy.sparse as sp
def joint_scores(method):
    RCC="kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad"
    PRO_scvi="prostate-cancer/Cleaned_Data/myeloid_integrated_final_label.h5ad"
    PRO_raw ="prostate-cancer/Cleaned_Data/myeloid_FINAL.h5ad"
    pairs=[("RCC",RCC),("prostate", PRO_scvi if method=="scvi" else PRO_raw)]
    ads=[]
    for tag,path in pairs:
        a=ad.read_h5ad(os.path.join(ROOT,path))
        if method=="scvi":
            M=a.layers["scvi_normalized"]; M=M.toarray() if sp.issparse(M) else np.asarray(M)
            sub=ad.AnnData(np.log1p(M),obs=a.obs.copy(),var=pd.DataFrame(index=a.var_names.astype(str)))
        else:
            sub=a.raw.to_adata() if (a.raw is not None and a.raw.n_vars>=a.n_vars) else a.copy()
            sub.obs=a.obs.copy()
        ads.append(sub)
    j=ad.concat(ads,join="inner")
    core_j=[g for g in CORE if g in j.var_names]
    dc.mt.aucell(j,pd.DataFrame({"source":"CORE","target":core_j,"weight":1.0}),tmin=3,verbose=False)
    print(f"  joint[{method}]: {j.n_vars} bg genes, CORE ({len(core_j)}/8)={core_j}")
    return pd.Series(j.obsm["score_aucell"]["CORE"].values,index=j.obs_names), core_j

def joint_pseudobulk(sj):
    pcj=pc.copy(); pcj["score"]=sj.reindex(pcj.index)
    q=(pcj.groupby(["Sample","cancer_type","patient_id","condition"],observed=True)
          .agg(score=("score","mean"),n_cells=("score","size")).reset_index())
    return q[q.n_cells>=MIN_CELLS].copy()
def joint_interaction(pbj):
    bj=pbj[pbj.condition.isin(["Tumor","Benign"])].copy()
    bj["condition"]=pd.Categorical(bj["condition"],["Benign","Tumor"]); bj["cancer_type"]=pd.Categorical(bj["cancer_type"],["prostate","RCC"])
    try: return mixedlm_coef(bj,form,"patient_id",coef)
    except Exception: return ols_coef(bj,form,coef,cluster="patient_id")

sj_scvi,core_scvi=joint_scores("scvi"); pbj=joint_pseudobulk(sj_scvi)
pbj.to_csv(os.path.join(TAB,"pseudobulk_joint_scVInorm.csv"),index=False)
e,l,h,p=joint_interaction(pbj); add("C2_interaction_scVInorm_ownBenign","RCC-skewed (integration arm: scVI-normalized, MAIN)",e,l,h,p,f"CORE={len(core_scvi)}/8")

sj_raw,core_raw=joint_scores("raw"); pbj_raw=joint_pseudobulk(sj_raw)
er,lr,hr,pr=joint_interaction(pbj_raw); add("C2_interaction_rawCommon_ownBenign","RCC-skewed (conservative common-space supplemental)",er,lr,hr,pr,f"CORE={len(core_raw)}/8")

def interaction_rcc_benign_variant(pbj, rcc_benign_df):
    pros=pbj[(pbj.cancer_type=="prostate")&pbj.condition.isin(["Tumor","Benign"])][["patient_id","cancer_type","condition","score"]]
    rcc_t=pbj[(pbj.cancer_type=="RCC")&(pbj.condition=="Tumor")][["patient_id","cancer_type","condition","score"]]
    rb=rcc_benign_df.copy(); rb["cancer_type"]="RCC"; rb["condition"]="Benign"
    df=pd.concat([pros,rcc_t,rb[["patient_id","cancer_type","condition","score"]]],ignore_index=True)
    df["condition"]=pd.Categorical(df["condition"],["Benign","Tumor"]); df["cancer_type"]=pd.Categorical(df["cancer_type"],["prostate","RCC"])
    return ols_coef(df,form,coef,cluster="patient_id")+(int((df.cancer_type=='RCC').sum()),)
bR=pbj[(pbj.cancer_type=="RCC")&(pbj.condition=="Benign")][["patient_id","score"]]
bP=pbj[(pbj.cancer_type=="prostate")&(pbj.condition=="Benign")][["patient_id","score"]]
bA=pd.concat([bR,bP]).groupby("patient_id",as_index=False)["score"].mean()
for nm,br in [("RCCbenign",bR),("PRObenign_forRCCarm",bP),("AVGbenign",bA)]:
    if len(br)==0: add(f"C2_benignRef_{nm}","RCC-skewed (benign-ref robustness)",np.nan,np.nan,np.nan,np.nan,"no benign rows pass gate"); continue
    e,l,h,p,n=interaction_rcc_benign_variant(pbj,br); add(f"C2_benignRef_{nm}","RCC-skewed (benign-ref robustness, scVI)",e,l,h,p,f"n_RCC_units={n}")

reg=b.copy()
dummies=pd.get_dummies(reg["dataset"],drop_first=True).astype(float)
reg["score"]=sm.OLS(reg["score"],sm.add_constant(dummies)).fit().resid
e,l,h,p=ols_coef(reg,form,coef,cluster="patient_id"); add("C2_interaction_datasetRegressed","RCC-skewed (integration arm3: dataset-regressed)",e,l,h,p)

comp=pc.assign(is_tt=pc.compartment.isin(["TAM","TIM"]))
tam_frac=(comp.groupby(["Sample","cancer_type","patient_id","condition"],observed=True)
              .agg(n_myeloid=("is_tt","size"),n_tt=("is_tt","sum")).reset_index())
tam_frac["tt_frac"]=tam_frac["n_tt"]/tam_frac["n_myeloid"]
tam_frac=tam_frac[tam_frac.n_myeloid>=MIN_CELLS]
tam_frac.to_csv(os.path.join(TAB,"tamtim_fraction_per_sample.csv"),index=False)

cf=tam_frac[tam_frac.condition.isin(["Tumor","Benign"])].copy()
cf["condition"]=pd.Categorical(cf["condition"],["Benign","Tumor"]); cf["cancer_type"]=pd.Categorical(cf["cancer_type"],["prostate","RCC"])
cf["asin"]=np.arcsin(np.sqrt(cf["tt_frac"].clip(0,1)))
e,l,h,p=ols_coef(cf,"asin ~ C(cancer_type, Treatment('prostate'))*C(condition, Treatment('Benign'))",coef,cluster="patient_id")
add("D1_composition_TAMTIMfrac_interaction","composition (TAM+TIM fraction, arcsin)",e,l,h,p,"is RCC-skew driven by more TAM+TIM?")

st=pc[pc.compartment.isin(["TAM","TIM"])].copy()
sst=(st.groupby(["Sample","cancer_type","patient_id","condition"],observed=True)
       .agg(score=("CLEC_LAM_aucell","mean"),n=("CLEC_LAM_aucell","size")).reset_index())
sst=sst[sst.n>=10]
tt_tum=sst[sst.condition=="Tumor"].copy(); tt_tum["cancer_type"]=pd.Categorical(tt_tum["cancer_type"],["prostate","RCC"])
e,l,h,p=ols_coef(tt_tum,"score ~ C(cancer_type, Treatment('prostate'))","C(cancer_type, Treatment('prostate'))[T.RCC]",cluster="patient_id")
add("D2_state_withinTAMTIM_tumor_RCCvsPro","state (within-TAM+TIM CORE, tumor-only)",e,l,h,p,
    f"RCC_tumor_n={(tt_tum.cancer_type=='RCC').sum()} pro_tumor_n={(tt_tum.cancer_type=='prostate').sum()}")

pst=sst[sst.cancer_type=="prostate"].copy()
if pst.condition.isin(["Benign"]).any():
    pst2=pst[pst.condition.isin(["Tumor","Benign"])].copy(); pst2["condition"]=pd.Categorical(pst2["condition"],["Benign","Tumor"])
    e,l,h,p=ols_coef(pst2,"score ~ C(condition, Treatment('Benign'))","C(condition, Treatment('Benign'))[T.Tumor]",cluster="patient_id")
    add("D2b_state_prostate_withinTAMTIM_TvsB","state (prostate within-TAM+TIM tumor vs benign)",e,l,h,p)

badj=b.merge(tam_frac[["Sample","tt_frac"]],on="Sample",how="left")
e,l,h,p=ols_coef(badj,form+" + tt_frac",coef,cluster="patient_id")
add("D3_interaction_adjTAMTIMfrac_RESIDUAL","residual (interaction | TAM+TIM frac) — GUARDED",e,l,h,p,
    "if attenuated: skew is compositional (=canalization); do NOT read as null")

res=pd.DataFrame(results)

conf=["C1_RCC_tumor_vs_benign_allmyeloid","C2_interaction_ownBenign_allmyeloid","C3b_fraction_interaction_arcsinLMM"]
res["fdr_confirmatory"]=np.nan; pv=res.loc[res.test.isin(conf)&res.p.notna(),"p"]
if len(pv): res.loc[pv.index,"fdr_confirmatory"]=multipletests(pv,method="fdr_bh")[1]
res.to_csv(os.path.join(TAB,"phase1_confirmatory_results.csv"),index=False)

print("=== TAM+TIM availability (structural) ==="); print(avail.to_string(index=False))
print("\n=== Phase 1 results ===")
with pd.option_context("display.width",220,"display.max_colwidth",64):
    print(res[["test","estimate","ci_low","ci_high","p","fdr_confirmatory","note"]].to_string(index=False))
print("\nLOO headline: est[%.4f,%.4f] max_p=%.3g (n_drop=%d)"%(loo.interaction_est.min(),loo.interaction_est.max(),loo.p.max(),len(loo)))
