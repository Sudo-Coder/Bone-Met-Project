#!/usr/bin/env python
"""41_braun_icb.py — Phase 4b: Braun/CheckMate ICB. Primary = nivolumab-treated OS/PFS/ORR/benefit;
exploratory = CM-025 nivo-vs-everolimus arm x module interaction. Control-contrast discipline.
Modules scored on Braun normalized expression (AUCell single-sample). Adjust for IMDC, Braun Myeloid
signature, Purity where available; report collinearity honestly. Run: envs/rcc_reinterp_venv/bin/python.
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, decoupler as dc
import statsmodels.formula.api as smf
from lifelines import CoxPHFitter
np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables")
B=os.path.join(ROOT,"resources/external_icb/braun_checkmate/processed")
MOD={"complement_C1Q":["C1QA","C1QB","C1QC"],"complement_C1Q_C3":["C1QA","C1QB","C1QC","C3"],
     "CLEC_LAM8":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"],
     "RCC_skew_CORE":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2"],
     "APOE_TREM2":["APOE","APOC1","TREM2","TYROBP"],"MERTK_GPNMB":["MERTK","GPNMB"],
     "panTAM":["CD68","CD163","MRC1","CSF1R","LYZ","AIF1","FCGR3A"],
     "Obradovic_TREM2":["TREM2","APOE","APOC1","C1QA","C1QB","C1QC","GPNMB","FOLR2","SPP1","CTSD","CD68"]}
expr=pd.read_csv(os.path.join(B,"expression_normalized.tsv"),sep="\t",index_col=0)
expr.index=[str(g).upper() for g in expr.index]
X=expr.T  # samples(RNA_ID) x genes
adT=ad.AnnData(X.values.astype(float),obs=pd.DataFrame(index=X.index),var=pd.DataFrame(index=X.columns))
net=pd.concat([pd.DataFrame({"source":k,"target":[g for g in v if g in adT.var_names],"weight":1.0}) for k,v in MOD.items()])
dc.mt.aucell(adT,net,tmin=2,verbose=False); SC=adT.obsm["score_aucell"]
sc=pd.DataFrame({k:SC[k].values for k in MOD if k in SC.columns},index=X.index); sc["RNA_ID"]=X.index
man=pd.read_csv(os.path.join(B,"sample_manifest.tsv"),sep="\t")
cov=pd.read_csv(os.path.join(B,"covariates.tsv"),sep="\t"); sv=pd.read_csv(os.path.join(B,"survival.tsv"),sep="\t"); rs=pd.read_csv(os.path.join(B,"response.tsv"),sep="\t")
D=sc.merge(man[["RNA_ID","SUBJID","Cohort","Arm"]],on="RNA_ID").merge(cov,on=["SUBJID","Cohort","Arm"],how="left").merge(sv,on="SUBJID",how="left").merge(rs,on="SUBJID",how="left")
D["OS_event"]=(D["OS_CNSR"]==0).astype(int); D["PFS_event"]=(D["PFS_CNSR"]==0).astype(int)
D["responder"]=D["ORR"].isin(["CR","PR","CRPR"]).astype(int)
D["benefit"]=D["Benefit"].map({"CB":1,"NCB":0}); D["IMDC_num"]=pd.to_numeric(D["IMDC"],errors="coerce")
for c in ["Myeloid","Purity","Age"]: D[c]=pd.to_numeric(D[c],errors="coerce")
def zc(x): x=pd.to_numeric(x,errors="coerce"); return (x-x.mean())/x.std()
print("Braun RNA samples:",len(D),"| arms:",dict(D.Arm.value_counts()),"| nivo cohorts:",dict(D[D.Arm=='NIVOLUMAB'].Cohort.value_counts()))

TEST=["complement_C1Q","complement_C1Q_C3","CLEC_LAM8","RCC_skew_CORE","APOE_TREM2","MERTK_GPNMB","panTAM","Obradovic_TREM2"]
nivo=D[D.Arm=="NIVOLUMAB"].copy()
rows=[]
def cox(df,ev,tm,m,adj):
    d=df.dropna(subset=[m,ev,tm]).copy(); d=d[d[tm]>0]
    if adj: d=d.dropna(subset=["Myeloid","Purity"])
    if len(d)<25 or d[ev].sum()<10: return dict(model="adj" if adj else "unadj",HR_per_SD=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d),events=int(d[ev].sum()) if len(d) else 0)
    d["z"]=zc(d[m]); cols=["z",tm,ev]
    if adj:
        d["z_my"]=zc(d["Myeloid"]); d["z_pu"]=zc(d["Purity"]); cols=["z","z_my","z_pu",tm,ev]
        if m=="panTAM": cols=["z","z_pu",tm,ev]
    try:
        cph=CoxPHFitter().fit(d[cols].rename(columns={tm:"T",ev:"E"}),"T","E"); r=cph.summary.loc["z"]
        return dict(model="adj(Myeloid+Purity)" if adj else "unadj",HR_per_SD=np.exp(r["coef"]),ci_low=np.exp(r["coef lower 95%"]),ci_high=np.exp(r["coef upper 95%"]),p=r["p"],n=len(d),events=int(d[ev].sum()))
    except Exception as e: return dict(model="adj" if adj else "unadj",HR_per_SD=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d),events=int(d[ev].sum()),note=str(e)[:30])
for m in TEST:
    for ep,ev,tm in [("OS","OS_event","OS"),("PFS","PFS_event","PFS")]:
        for adj in (False,True):
            r=cox(nivo,ev,tm,m,adj); r.update(dict(module=m,endpoint=ep)); rows.append(r)
# response logistic
for m in TEST:
    for oc in ["responder","benefit"]:
        d=nivo.dropna(subset=[m,oc]).copy()
        if len(d)<25 or d[oc].sum()<5: rows.append(dict(module=m,endpoint=oc,model="logit",HR_per_SD=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d),events=int(d[oc].sum()) if len(d) else 0)); continue
        d["z"]=zc(d[m])
        try:
            fit=smf.logit(f"{oc} ~ z",d).fit(disp=0); ci=fit.conf_int().loc["z"]
            rows.append(dict(module=m,endpoint=oc,model="logit(OR/SD)",HR_per_SD=np.exp(fit.params["z"]),ci_low=np.exp(ci[0]),ci_high=np.exp(ci[1]),p=fit.pvalues["z"],n=len(d),events=int(d[oc].sum())))
        except Exception as e: rows.append(dict(module=m,endpoint=oc,model="logit",HR_per_SD=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d),events=int(d[oc].sum())))
R=pd.DataFrame(rows); R.to_csv(os.path.join(TAB,"phase4_braun_nivo.csv"),index=False)

# CM-025 nivo-vs-evero arm interaction (exploratory)
cm=D[D.Cohort=="CM-025"].copy(); cm["arm_nivo"]=(cm.Arm=="NIVOLUMAB").astype(int)
h=[]
for m in ["complement_C1Q","CLEC_LAM8","RCC_skew_CORE"]:
    d=cm.dropna(subset=[m,"OS","OS_event"]).copy(); d=d[d.OS>0]; d["z"]=zc(d[m]); d["zxarm"]=d["z"]*d["arm_nivo"]
    try:
        cph=CoxPHFitter().fit(d[["z","arm_nivo","zxarm","OS","OS_event"]].rename(columns={"OS":"T","OS_event":"E"}),"T","E")
        r=cph.summary.loc["zxarm"]; h.append(dict(module=m,term="module x nivo-arm (OS)",HR=np.exp(r["coef"]),ci_low=np.exp(r["coef lower 95%"]),ci_high=np.exp(r["coef upper 95%"]),p=r["p"],n=len(d),n_nivo=int(d.arm_nivo.sum()),n_evero=int((1-d.arm_nivo).sum())))
    except Exception as e: h.append(dict(module=m,term="interaction",HR=np.nan,p=np.nan,note=str(e)[:40]))
H=pd.DataFrame(h); H.to_csv(os.path.join(TAB,"phase4_braun_cm025_armINT.csv"),index=False)
pd.set_option("display.width",230)
print("\n=== Braun nivolumab-treated (primary) ===")
print(R[["module","endpoint","model","HR_per_SD","ci_low","ci_high","p","n","events"]].round(3).to_string(index=False))
print("\n=== CM-025 nivo-vs-evero arm x module interaction (EXPLORATORY, underpowered) ===")
print(H.round(3).to_string(index=False))
