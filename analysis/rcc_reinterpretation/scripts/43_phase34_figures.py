#!/usr/bin/env python
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
FIG=os.path.join(ROOT,"analysis/rcc_reinterpretation/outputs/figures"); TAB=os.path.join(ROOT,"analysis/rcc_reinterpretation/outputs/tables")
tcga=pd.read_csv(os.path.join(TAB,"phase4_tcga_kirc_cox.csv")); p3=pd.read_csv(os.path.join(TAB,"phase3_adjusted_associations.csv"))
MODS=["complement_C1Q","RCC_skew_CORE","CLEC_LAM8","APOE_TREM2","MERTK_GPNMB","panTAM"]

d=tcga[(tcga.endpoint=="OS")].copy()
d=d[d.model.str.startswith("adj") | ((d.module.isin(["panTAM"]))&(d.model=="unadj"))]
d=d[d.module.isin(MODS)].drop_duplicates("module")
d=d.set_index("module").reindex(MODS).reset_index()
fig,ax=plt.subplots(figsize=(7,4)); y=np.arange(len(d))[::-1]
ax.errorbar(d.HR_per_SD,y,xerr=[d.HR_per_SD-d.ci_low,d.ci_high-d.HR_per_SD],fmt="o",color="#333",capsize=3)
ax.axvline(1,color="r",ls="--",lw=1); ax.set_yticks(y); ax.set_yticklabels([f"{m}  (p={p:.3f})" for m,p in zip(d.module,d.p)],fontsize=8)
ax.set_xlabel("TCGA-KIRC OS  HR per SD (adjusted: panTAM+Obradovic+immune)"); ax.set_title("Phase 4 — TCGA-KIRC OS (continuous Cox, adjusted)",fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"figP34_1_tcga_os_forest.png"),dpi=300); plt.close()

outs=["CD8_exhaustion","cytotoxicity","Treg_fraction","MHC_II_APC"]
sub=p3[p3.module.isin(["complement_C1Q","panTAM","SPP1_TAM","MERTK_GPNMB"])&p3.outcome.isin(outs)]
piv=sub.pivot_table(index="outcome",columns="module",values="coef").reindex(outs)
fig,ax=plt.subplots(figsize=(8,4)); x=np.arange(len(outs)); w=0.2
cols={"complement_C1Q":"#2E7D32","panTAM":"#9E9E9E","SPP1_TAM":"#C99B38","MERTK_GPNMB":"#7B8794"}
for i,(m,c) in enumerate(cols.items()):
    if m in piv.columns: ax.bar(x+(i-1.5)*w,piv[m],w,label=m,color=c)
ax.axhline(0,color="k",lw=.6); ax.set_xticks(x); ax.set_xticklabels(outs,fontsize=8,rotation=20)
ax.set_ylabel("adjusted std. coef"); ax.legend(fontsize=7); ax.set_title("Phase 3 — immune-evasion contrast (complement vs generic/other TAM)",fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"figP34_2_immune_evasion_contrast.png"),dpi=300); plt.close()
print("wrote figP34_1, figP34_2 (300dpi)")
