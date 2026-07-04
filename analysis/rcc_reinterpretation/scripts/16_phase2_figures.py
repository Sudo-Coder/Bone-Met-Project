#!/usr/bin/env python
"""16_phase2_figures.py — Phase 2 figures (300 dpi). Run: envs/rcc_reinterp_venv/bin/python."""
import os, warnings, ast
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
FIG=os.path.join(ROOT,"analysis/rcc_reinterpretation/outputs/figures"); os.makedirs(FIG,exist_ok=True)
TAB=os.path.join(ROOT,"analysis/rcc_reinterpretation/outputs/tables")
ax=pd.read_csv(os.path.join(TAB,"axis_support_table.csv"))
lr=pd.read_csv(os.path.join(TAB,"axis_lr_expression.csv"))
cc=pd.read_csv(os.path.join(TAB,"cellchat_tam_LR_tumor.csv"))
nn=pd.read_csv(os.path.join(TAB,"nichenet_prespecified_ranks.csv"))
VC={"SUPPORTED_TUMOR_SPECIFIC":"#2E7D32","AMBIGUOUS":"#C99B38","PRESENT_NOT_TUMOR_GAINED":"#7B8794",
    "NOT_TUMOR_GAINED":"#B0553F","UNSUPPORTED":"#9E9E9E"}

# Fig P2-1: axis support summary
fig,axs=plt.subplots(1,3,figsize=(13,4.2),gridspec_kw={"width_ratios":[1.1,1,1]})
a=ax.copy(); a["short"]=a["axis"].str.replace("complement_","C").str.replace("_TREM2","→TREM2").str.replace("_MERTK","→MERTK").str.replace("_TGFBR","→TGFBR")
y=np.arange(len(a))[::-1]
axs[0].barh(y,a["lig_delta_TvsB"],color=["#2E7D32" if v else "#B0553F" for v in a["A_tumor_elevated"]])
axs[0].axvline(0,color="k",lw=.6); axs[0].set_yticks(y); axs[0].set_yticklabels(a["short"],fontsize=8)
axs[0].set_title("(A) ligand Δ tumor−benign\n(green=tumor-elevated,p<0.05)",fontsize=8); axs[0].set_xlabel("Δ log-norm")
axs[1].barh(y,a["B_cellchat_maxprob"].fillna(0),color="#4C72B0")
axs[1].set_yticks(y); axs[1].set_yticklabels([]); axs[1].set_title("(B) CellChat max prob\nsender→TAM",fontsize=8); axs[1].set_xlabel("comm. prob")
rank=a["C_nichenet_rank"].fillna(a["C_nichenet_rank"].max()+10)
axs[2].barh(y,rank,color=["#2E7D32" if (v is True or v=="True") else "#7B8794" for v in a["C_top20"]])
axs[2].set_yticks(y); axs[2].set_yticklabels([]); axs[2].invert_xaxis(); axs[2].set_title("(C) NicheNet rank\n(of 116; lower=better)",fontsize=8); axs[2].set_xlabel("rank")
for i,v in zip(y,a["verdict"]): axs[2].text(rank.iloc[len(y)-1-i]*0+2,i,"",fontsize=6)
# verdict legend
handles=[plt.Rectangle((0,0),1,1,color=c) for c in VC.values()]
fig.legend(handles,VC.keys(),loc="lower center",ncol=5,fontsize=6.5,frameon=False,bbox_to_anchor=(0.5,-0.03))
plt.suptitle("Phase 2 — pre-specified axis support (predicted signaling axes)",fontsize=11)
plt.tight_layout(rect=[0,0.04,1,0.96]); plt.savefig(os.path.join(FIG,"figP2_1_axis_support.png"),dpi=300); plt.close()

# Fig P2-2: CellChat top sender->TAM_CLEC_LAM dotplot
cl=cc[cc.target=="TAM_CLEC_LAM"].sort_values("prob",ascending=False).head(18).copy()
cl["pair"]=cl["source"]+" → "+cl["ligand"]+"→"+cl["receptor"]
fig,axd=plt.subplots(figsize=(7.5,6))
yy=np.arange(len(cl))[::-1]
sc=axd.scatter(cl["prob"],yy,s=120,c=cl["prob"],cmap="viridis")
axd.set_yticks(yy); axd.set_yticklabels(cl["pair"],fontsize=7.5); axd.set_xlabel("communication probability")
axd.set_title("CellChat: top predicted axes → TAM_CLEC_LAM (tumor niche)",fontsize=10)
plt.colorbar(sc,label="prob"); plt.tight_layout(); plt.savefig(os.path.join(FIG,"figP2_2_cellchat_clec_lam.png"),dpi=300); plt.close()

# Fig P2-3: NicheNet pre-specified ranks
nn2=nn.sort_values("rank"); fig,axn=plt.subplots(figsize=(6,3.4))
axn.barh(np.arange(len(nn2))[::-1], nn2["aupr_corrected"],
         color=["#2E7D32" if r<=20 else "#7B8794" for r in nn2["rank"]])
axn.set_yticks(np.arange(len(nn2))[::-1]); axn.set_yticklabels([f"{t} (r{int(r)})" for t,r in zip(nn2.test_ligand,nn2["rank"])],fontsize=8)
axn.set_xlabel("NicheNet aupr_corrected"); axn.set_title("NicheNet ligand→CLEC_LAM-program activity\n(green=top20 of 116)",fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"figP2_3_nichenet_ranks.png"),dpi=300); plt.close()

# Fig P2-4: ligand tumor vs benign
fig,axl=plt.subplots(figsize=(7,3.6)); x=np.arange(len(lr)); w=0.38
axl.bar(x-w/2,lr["ligand_benign_mean"],w,label="benign senders",color="#4C72B0")
axl.bar(x+w/2,lr["ligand_tumor_mean"],w,label="tumor senders",color="#C44E52")
axl.set_xticks(x); axl.set_xticklabels(lr["ligand"],rotation=45,ha="right",fontsize=8)
axl.set_ylabel("mean log-norm (patient)"); axl.legend(fontsize=8); axl.set_title("Pre-specified ligand expression in senders (tumor vs benign)",fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"figP2_4_ligand_tumor_vs_benign.png"),dpi=300); plt.close()
print("wrote figP2_1..4 (300 dpi)")
