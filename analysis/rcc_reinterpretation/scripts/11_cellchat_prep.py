#!/usr/bin/env python
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py, anndata as ad, decoupler as dc
import scipy.sparse as sp, scipy.io as sio
np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables")
OUTD=os.path.join(TAB,"cellchat_input"); os.makedirs(OUTD,exist_ok=True)
WS=os.path.join(ROOT,"kidney-cancer/Cleaned_Data/integrated-with-stromal.h5ad")
CORE=["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"]

f=h5py.File(WS,"r")
def idx(g): k=g.attrs.get("_index","_index"); return k.decode() if isinstance(k,bytes) else k
names=np.array([x.decode() for x in f["obs"][idx(f["obs"])][:]])
vnames=np.array([x.decode() for x in f["var"][idx(f["var"])][:]])
Xg=f["X"]; X=sp.csr_matrix((Xg["data"][:],Xg["indices"][:],Xg["indptr"][:]),shape=tuple(Xg.attrs["shape"]))
f.close()

lab=pd.read_csv(os.path.join(TAB,"withstromal_labels.csv")).set_index("cell").reindex(names)

A=ad.AnnData(X.copy(),obs=pd.DataFrame(index=names),var=pd.DataFrame(index=vnames))
core=[g for g in CORE if g in set(vnames)]
dc.mt.aucell(A,pd.DataFrame({"source":"CORE","target":core,"weight":1.0}),tmin=3,verbose=False)
sc_core=A.obsm["score_aucell"]["CORE"].values
tt=lab["compartment"].isin(["TAM","TIM"]).values
cut=np.quantile(sc_core[tt],0.80) if tt.sum() else np.inf
group=lab["compartment"].astype(str).values.copy()
is_tam=lab["compartment"].values=="TAM"
group[is_tam & (sc_core>=cut)]="TAM_CLEC_LAM"
group[is_tam & (sc_core<cut)]="TAM_other"
lab["cc_group"]=group; lab["CORE_aucell"]=sc_core

KEEP=["Tumor","MSC","Pericyte","Endothelial","Osteoclast","TAM_CLEC_LAM","TAM_other","TIM","T_NK"]
def export(niche, cond_mask):
    idxs=np.where(cond_mask & lab["cc_group"].isin(KEEP).values)[0]
    d=os.path.join(OUTD,niche); os.makedirs(d,exist_ok=True)
    sub=X[idxs].T.tocsc()
    sio.mmwrite(os.path.join(d,"expr.mtx"), sub)
    pd.Series(vnames).to_csv(os.path.join(d,"genes.txt"),index=False,header=False)
    pd.Series(names[idxs]).to_csv(os.path.join(d,"cells.txt"),index=False,header=False)
    m=lab.iloc[idxs][["Sample","condition","cc_group","CORE_aucell","transfer_conf"]].copy()
    m.index=names[idxs]; m.to_csv(os.path.join(d,"meta.csv"))
    print(f"[{niche}] {len(idxs)} cells | groups:", dict(m["cc_group"].value_counts()))
    return m

tumor_mask=lab["condition"].values=="Tumor"
benign_mask=np.isin(lab["condition"].values,["Benign_stroma","Benign_immune"])
mt=export("tumor",tumor_mask)
mb=export("benign",benign_mask)
print("\nbenign receiver check — TAM/TIM in benign:",
      int(mb["cc_group"].isin(["TAM_CLEC_LAM","TAM_other","TIM"]).sum()),
      "(=> benign->TAM CellChat unpopulated; tumor-specificity from LR expression)")
print("wrote cellchat_input/{tumor,benign}/")
