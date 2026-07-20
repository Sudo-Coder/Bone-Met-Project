#!/usr/bin/env python
import os, warnings, gzip
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, decoupler as dc
from lifelines import CoxPHFitter
np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables")
EXP=os.path.join(ROOT,"resources/tcga/KIRC_HiSeqV2.gz"); SURV=os.path.join(ROOT,"resources/tcga/KIRC_survival.txt")
MOD={"complement_C1Q":["C1QA","C1QB","C1QC"],"complement_C1Q_C3":["C1QA","C1QB","C1QC","C3"],
     "CLEC_LAM8":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"],
     "RCC_skew_CORE":["C1QA","C1QB","C1QC","APOE","APOC1","TREM2"],
     "APOE_TREM2":["APOE","APOC1","TREM2","TYROBP"],"MERTK_GPNMB":["MERTK","GPNMB"],
     "panTAM":["CD68","CD163","MRC1","CSF1R","LYZ","AIF1","FCGR3A"],
     "Obradovic_TREM2":["TREM2","APOE","APOC1","C1QA","C1QB","C1QC","GPNMB","FOLR2","SPP1","CTSD","CD68"],
     "total_immune":["PTPRC","CD3D","CD8A","CD4","CD19","NKG7","CD68"]}

expr=pd.read_csv(EXP,sep="\t",index_col=0)
expr.index=[str(g).upper() for g in expr.index]
samples=[c for c in expr.columns if c[-2:] in ("01","05")]
X=expr[samples].T
adT=ad.AnnData(X.values.astype(float),obs=pd.DataFrame(index=X.index),var=pd.DataFrame(index=X.columns))
net=pd.concat([pd.DataFrame({"source":k,"target":[g for g in v if g in adT.var_names],"weight":1.0}) for k,v in MOD.items()])
dc.mt.aucell(adT,net,tmin=2,verbose=False); SC=adT.obsm["score_aucell"]
sc=pd.DataFrame({k:SC[k].values for k in MOD if k in SC.columns},index=X.index)

sv=pd.read_csv(SURV,sep="\t")
sv["bcr"]=sv["sample"].str[:15]
sc["bcr"]=[s[:15] for s in sc.index]
D=sc.merge(sv,on="bcr",how="inner")
print("TCGA tumor samples scored:",len(sc),"| merged w/ survival:",len(D))
def zc(x): x=pd.to_numeric(x,errors="coerce"); return (x-x.mean())/x.std()

endpoints={"OS":("OS","OS.time"),"PFI":("PFI","PFI.time"),"DFI":("DFI","DFI.time")}
TEST=["complement_C1Q","complement_C1Q_C3","CLEC_LAM8","RCC_skew_CORE","APOE_TREM2","MERTK_GPNMB","panTAM","Obradovic_TREM2"]
rows=[]
for ep,(ev,tm) in endpoints.items():
    d0=D.dropna(subset=[ev,tm]).copy(); d0=d0[d0[tm]>0]
    for m in TEST:
        d=d0.dropna(subset=[m]).copy();
        if len(d)<30 or d[ev].sum()<10: rows.append(dict(endpoint=ep,module=m,model="unadj",HR_per_SD=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d),events=int(d[ev].sum()))); continue
        d["z"]=zc(d[m])

        try:
            cph=CoxPHFitter().fit(d[["z",tm,ev]].rename(columns={tm:"T",ev:"E"}),"T","E")
            r=cph.summary.loc["z"]; rows.append(dict(endpoint=ep,module=m,model="unadj",HR_per_SD=np.exp(r["coef"]),ci_low=np.exp(r["coef lower 95%"]),ci_high=np.exp(r["coef upper 95%"]),p=r["p"],n=len(d),events=int(d[ev].sum())))
        except Exception as e: rows.append(dict(endpoint=ep,module=m,model="unadj",HR_per_SD=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d),events=int(d[ev].sum())))

        try:
            d["z_pan"]=zc(d["panTAM"]); d["z_obr"]=zc(d["Obradovic_TREM2"]); d["z_imm"]=zc(d["total_immune"])
            cols=["z","z_pan","z_obr","z_imm",tm,ev]
            if m in ("panTAM","Obradovic_TREM2"): cols=["z","z_imm",tm,ev]
            cph=CoxPHFitter().fit(d[cols].rename(columns={tm:"T",ev:"E"}),"T","E")
            r=cph.summary.loc["z"]; rows.append(dict(endpoint=ep,module=m,model="adj(panTAM+Obradovic+immune)",HR_per_SD=np.exp(r["coef"]),ci_low=np.exp(r["coef lower 95%"]),ci_high=np.exp(r["coef upper 95%"]),p=r["p"],n=len(d),events=int(d[ev].sum())))
        except Exception as e: rows.append(dict(endpoint=ep,module=m,model="adj",HR_per_SD=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,n=len(d),events=int(d[ev].sum()),note=str(e)[:40]))
R=pd.DataFrame(rows); R.to_csv(os.path.join(TAB,"phase4_tcga_kirc_cox.csv"),index=False)
pd.set_option("display.width",220)
print("\n=== TCGA-KIRC continuous Cox (HR per SD) ===")
print(R[["endpoint","module","model","HR_per_SD","ci_low","ci_high","p","n","events"]].round(3).to_string(index=False))
print("\nUnavailable covariates in staged TCGA file: stage, grade, age, sex, tumor purity (parsimonious models; noted).")
