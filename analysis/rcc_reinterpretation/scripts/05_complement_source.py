#!/usr/bin/env python
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables")
WS=os.path.join(ROOT,"kidney-cancer/Cleaned_Data/integrated-with-stromal.h5ad")
GENES=["C1QA","C1QB","C1QC","C3","C1R","C1S","CFB","SERPING1"]
CLASSES=["Tumor","TAM_CLEC_LAM","TAM_other","TIM","MSC","Pericyte","Endothelial","Osteoclast"]

lab=pd.read_csv(os.path.join(TAB,"withstromal_labels.csv")).set_index("cell")
ccm=pd.concat([pd.read_csv(os.path.join(TAB,"cellchat_input",n,"meta.csv"),index_col=0) for n in ["tumor","benign"]])
lab["cc_group"]=ccm["cc_group"].reindex(lab.index); lab["cc_group"]=lab["cc_group"].fillna(lab["compartment"])
def pid(s):
    s=str(s).split("_",1)[-1]
    for k in ["-Tumor",".","-Benign","-Involve","-Noninvolved"]: s=s.split(k)[0]
    return s.replace("RCC-","")
lab["pid"]=lab["Sample"].map(pid)

f=h5py.File(WS,"r")
def idx(g): k=g.attrs.get("_index","_index"); return k.decode() if isinstance(k,bytes) else k
names=np.array([x.decode() for x in f["obs"][idx(f["obs"])][:]])
vnames=list(np.array([x.decode() for x in f["var"][idx(f["var"])][:]])); gi={g:i for i,g in enumerate(vnames)}
Xg=f["X"]; X=sp.csr_matrix((Xg["data"][:],Xg["indices"][:],Xg["indptr"][:]),shape=tuple(Xg.attrs["shape"])).tocsc()
f.close()
lab=lab.reindex(names)
tumor=lab["condition"].values=="Tumor"

rows=[]
for g in GENES:
    if g not in gi:
        for cl in CLASSES: rows.append(dict(gene=g,compartment=cl,mean_lognorm=np.nan,frac_expr=np.nan,n=0))
        continue
    v=np.asarray(X[:,gi[g]].todense()).ravel()
    for cl in CLASSES:
        m=(lab["cc_group"].values==cl)&tumor
        if m.sum()==0: rows.append(dict(gene=g,compartment=cl,mean_lognorm=np.nan,frac_expr=np.nan,n=0)); continue

        df=pd.DataFrame({"x":v[m],"pid":lab["pid"].values[m]})
        pm=df.groupby("pid")["x"].mean()
        rows.append(dict(gene=g,compartment=cl,mean_lognorm=float(pm.mean()),
                         frac_expr=float((v[m]>0).mean()),n=int(m.sum())))
src=pd.DataFrame(rows); src.to_csv(os.path.join(TAB,"complement_source_by_class.csv"),index=False)
piv=src.pivot_table(index="compartment",columns="gene",values="mean_lognorm").reindex(CLASSES).dropna(axis=1,how="all")
print("=== C1q/C3 mean log-norm by sender class (tumor niche, patient-level) ===")
print(piv.round(2).to_string())
print("\n=== fraction expressing ===")
pf0=src.pivot_table(index="compartment",columns="gene",values="frac_expr").reindex(CLASSES); pf=pf0[[c for c in ["C1QA","C1QB","C1QC","C3"] if c in pf0.columns]]
print(pf.round(2).to_string())

for g in ["C1QA","C1QB","C1QC","C3"]:
    sub=src[(src.gene==g)&src.compartment.isin(CLASSES)].dropna(subset=["mean_lognorm"])
    if len(sub):
        top=sub.sort_values("mean_lognorm",ascending=False).iloc[0]
        print(f"{g}: top producer = {top.compartment} (mean {top.mean_lognorm:.2f})")
