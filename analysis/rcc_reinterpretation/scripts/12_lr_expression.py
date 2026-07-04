#!/usr/bin/env python
"""12_lr_expression.py — Phase 2: patient-level ligand/receptor expression for the pre-specified axes.

Because benign has ~0 TAM receivers (TAMs tumor-restricted), tumor-specificity of sender->TAM axes is
derived HERE (not from benign CellChat): for each axis, ligand expression in its SENDER compartment
(tumor vs benign, patient-level Δ + Mann-Whitney) and receptor expression/detection in TAM (tumor).
Wording downstream: "predicted axes operating in the tumor niche"; ligand tumor-elevation = tumor-specificity.

Run: envs/rcc_reinterp_venv/bin/python. Seed 0.
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables")
WS=os.path.join(ROOT,"kidney-cancer/Cleaned_Data/integrated-with-stromal.h5ad")

# axis: (name, ligand, ligand_senders, receptor(s) on TAM, note)
AXES=[
 ("complement_C1QA","C1QA",["TAM_CLEC_LAM","TAM_other","MSC"],["C1QBP"],"C1q macrophage-produced; autocrine/paracrine"),
 ("complement_C1QB","C1QB",["TAM_CLEC_LAM","TAM_other"],["C1QBP"],"macrophage-produced"),
 ("complement_C1QC","C1QC",["TAM_CLEC_LAM","TAM_other"],["C1QBP"],"macrophage-produced"),
 ("complement_C3","C3",["MSC","Pericyte","Tumor"],["C3AR1","ITGAM","ITGB2"],"stromal/tumor C3 -> TAM CR"),
 ("APOE_TREM2","APOE",["TAM_CLEC_LAM","TAM_other","Tumor"],["TREM2","LRP1","SORL1"],"APOE macrophage+tumor derived; both directions"),
 ("GAS6_MERTK","GAS6",["MSC","Pericyte","Endothelial","TAM_other"],["MERTK","AXL"],"efferocytosis"),
 ("PROS1_MERTK","PROS1",["Tumor","MSC","Endothelial"],["MERTK","AXL"],"efferocytosis"),
 ("TGFB1_TGFBR","TGFB1",["Tumor","MSC","Pericyte","TAM_other"],["TGFBR1","TGFBR2"],"APOE/TGFb/ICB-resistance link"),
]
LIG=sorted({a[1] for a in AXES}); REC=sorted({r for a in AXES for r in a[3]})
GENES=sorted(set(LIG+REC))

lab=pd.read_csv(os.path.join(TAB,"withstromal_labels.csv"))
# re-derive cc_group with CLEC_LAM split (mirror 11): need CORE score — reuse from cellchat meta
ccm=pd.concat([pd.read_csv(os.path.join(TAB,"cellchat_input",n,"meta.csv"),index_col=0) for n in ["tumor","benign"]])
lab=lab.set_index("cell"); lab["cc_group"]=ccm["cc_group"].reindex(lab.index); lab["cc_group"]=lab["cc_group"].fillna(lab["compartment"])

f=h5py.File(WS,"r")
def idx(g): k=g.attrs.get("_index","_index"); return k.decode() if isinstance(k,bytes) else k
names=np.array([x.decode() for x in f["obs"][idx(f["obs"])][:]])
vnames=list(np.array([x.decode() for x in f["var"][idx(f["var"])][:]]))
gi={g:i for i,g in enumerate(vnames)}
Xg=f["X"]; X=sp.csr_matrix((Xg["data"][:],Xg["indices"][:],Xg["indptr"][:]),shape=tuple(Xg.attrs["shape"])).tocsc()
f.close()
lab=lab.reindex(names)
def gene_vec(g): return np.asarray(X[:,gi[g]].todense()).ravel() if g in gi else None
def patient_id(s):
    s=str(s).split("_",1)[-1]
    for kw in ["-Tumor",".","-Benign","-Involve","-Noninvolved"]: s=s.split(kw)[0]
    return s.replace("RCC-","")
lab["pid"]=lab["Sample"].map(patient_id)
lab["is_tumor"]=lab["condition"].eq("Tumor")
lab["is_benign"]=lab["condition"].isin(["Benign_stroma","Benign_immune"])

rows=[]
for name,lig,senders,recs,note in AXES:
    lv=gene_vec(lig)
    # ligand tumor vs benign in its sender compartments (patient-level mean)
    m=lab["cc_group"].isin(senders).values & (lv is not None)
    if lv is not None:
        df=pd.DataFrame({"x":lv,"pid":lab["pid"].values,"cond":np.where(lab["is_tumor"],"Tumor",np.where(lab["is_benign"],"Benign","other"))})[m]
        pm=df[df.cond.isin(["Tumor","Benign"])].groupby(["pid","cond"],observed=True)["x"].mean().reset_index()
        t=pm[pm.cond=="Tumor"]["x"]; b=pm[pm.cond=="Benign"]["x"]
        dlt=float(t.mean()-b.mean()) if len(t) and len(b) else np.nan
        try: _,pv=mannwhitneyu(t,b) if len(t) and len(b) else (np.nan,np.nan)
        except Exception: pv=np.nan
    else: dlt,pv,t,b=np.nan,np.nan,[],[]
    # receptor detection in tumor TAM (CLEC_LAM + other)
    tam=lab["cc_group"].isin(["TAM_CLEC_LAM","TAM_other"]).values & lab["is_tumor"].values
    recinfo={}
    for r in recs:
        rv=gene_vec(r)
        recinfo[r]=round(float((rv[tam]>0).mean()),3) if rv is not None else None
    rows.append(dict(axis=name,ligand=lig,ligand_senders=";".join(senders),
                     ligand_tumor_mean=round(float(np.mean(t)),3) if len(t) else np.nan,
                     ligand_benign_mean=round(float(np.mean(b)),3) if len(b) else np.nan,
                     ligand_delta_TvsB=round(dlt,3) if dlt==dlt else np.nan, ligand_p=pv,
                     receptors=";".join(recs), receptor_detect_in_tumorTAM=str(recinfo), note=note))
res=pd.DataFrame(rows); res.to_csv(os.path.join(TAB,"axis_lr_expression.csv"),index=False)
pd.set_option("display.width",240,"display.max_colwidth",40)
print(res[["axis","ligand","ligand_tumor_mean","ligand_benign_mean","ligand_delta_TvsB","ligand_p","receptor_detect_in_tumorTAM"]].to_string(index=False))
print("\nwrote axis_lr_expression.csv")
