#!/usr/bin/env python
import os, warnings, json
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

np.random.seed(0)
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
FIG=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","figures"); os.makedirs(FIG,exist_ok=True)
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables")
CORE=["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"]
COND_ORDER=["Benign","Distal","Involved","Tumor"]; COND_COL={"Benign":"#4C72B0","Distal":"#55A868","Involved":"#C99B38","Tumor":"#C44E52"}
pc=pd.read_parquet(os.path.join(TAB,"clec_lam_percell.parquet"))
res=pd.read_csv(os.path.join(TAB,"phase1_confirmatory_results.csv"))
OBJS={"RCC":"kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad",
      "prostate":"prostate-cancer/Cleaned_Data/myeloid_FINAL.h5ad"}

fig,axes=plt.subplots(2,2,figsize=(11,10))
for i,(tag,p) in enumerate(OBJS.items()):
    a=ad.read_h5ad(os.path.join(ROOT,p)); um=a.obsm["X_umap"]
    s=pc.loc[a.obs_names,"CLEC_LAM_aucell"].values
    cond=pc.loc[a.obs_names,"condition"].values
    sca=axes[0,i].scatter(um[:,0],um[:,1],c=s,s=3,cmap="magma",vmax=np.quantile(s,0.99))
    axes[0,i].set_title(f"{tag}: CLEC_LAM CORE (AUCell)"); plt.colorbar(sca,ax=axes[0,i],shrink=.7)
    for c in COND_ORDER:
        m=cond==c; axes[1,i].scatter(um[m,0],um[m,1],s=3,c=COND_COL[c],label=f"{c} ({m.sum()})")
    axes[1,i].legend(markerscale=3,fontsize=7); axes[1,i].set_title(f"{tag}: condition")
    for ax in (axes[0,i],axes[1,i]): ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); plt.savefig(os.path.join(FIG,"fig1_umap_core_score.png"),dpi=300); plt.close()

rep=json.load(open(os.path.join(TAB,"phase1_scoring_report.json")))
rows=[]
for tag in OBJS:
    for grp in ["TAM","TIM"]:
        d=rep["core_detection"][tag][grp]
        for g in CORE: rows.append(dict(cancer=tag,comp=grp,gene=g,det=d.get(g,np.nan)))
dd=pd.DataFrame(rows); piv=dd.pivot_table(index=["cancer","comp"],columns="gene",values="det")[CORE]
fig,ax=plt.subplots(figsize=(8,3.2)); im=ax.imshow(piv.values,cmap="viridis",vmin=0,vmax=1,aspect="auto")
ax.set_xticks(range(len(CORE))); ax.set_xticklabels(CORE,rotation=45,ha="right")
ax.set_yticks(range(len(piv))); ax.set_yticklabels([f"{a}|{b}" for a,b in piv.index])
for i in range(len(piv)):
    for j in range(len(CORE)): ax.text(j,i,f"{piv.values[i,j]:.2f}",ha="center",va="center",color="w",fontsize=7)
plt.colorbar(im,label="fraction expressing"); ax.set_title("CORE gene detection")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"fig2_core_detection.png"),dpi=300); plt.close()

pb=pd.read_csv(os.path.join(TAB,"pseudobulk_all_myeloid.csv")); pb=pb[pb.ok]
fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
for i,tag in enumerate(OBJS):
    sub=pb[pb.cancer_type==tag]; data=[sub[sub.condition==c]["score"].values for c in COND_ORDER]
    bp=axes[i].boxplot(data,labels=COND_ORDER,patch_artist=True,showfliers=False)
    for patch,c in zip(bp["boxes"],COND_ORDER): patch.set_facecolor(COND_COL[c]); patch.set_alpha(.6)
    for j,c in enumerate(COND_ORDER):
        y=sub[sub.condition==c]["score"].values; axes[i].scatter(np.full(len(y),j+1)+np.random.uniform(-.1,.1,len(y)),y,s=18,c="k",zorder=3)
    axes[i].set_title(f"{tag} — all-myeloid CORE score / patient"); axes[i].set_ylabel("mean AUCell CORE")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"fig3_patient_boxplots.png"),dpi=300); plt.close()

show=res[res.estimate.notna() & res.ci_low.notna() & (res.test.str.startswith(("C1","C2","C3b")))].copy()
show=show[~show.test.str.contains("binomial")]
fig,ax=plt.subplots(figsize=(9,0.5*len(show)+1))
y=np.arange(len(show))[::-1]
ax.errorbar(show.estimate,y,xerr=[show.estimate-show.ci_low,show.ci_high-show.estimate],fmt="o",color="#333",capsize=3)
ax.axvline(0,color="r",ls="--",lw=1); ax.set_yticks(y); ax.set_yticklabels(show.test,fontsize=8)
ax.set_xlabel("estimate (score units / arcsin) with 95% CI"); ax.set_title("Phase 1 confirmatory + robustness")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"fig4_forest.png"),dpi=300); plt.close()

gene_rows=[]
for tag,p in OBJS.items():
    a=ad.read_h5ad(os.path.join(ROOT,p))
    s=a.raw.to_adata() if (a.raw is not None and a.raw.n_vars>=a.n_vars) else a
    obs=pc.loc[a.obs_names]
    import scipy.sparse as sp
    for g in CORE:
        if g not in s.var_names: gene_rows.append(dict(cancer=tag,gene=g,delta=np.nan)); continue
        x=s[:,g].X; x=x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
        df=pd.DataFrame({"x":x,"cond":obs["condition"].values,"pid":obs["patient_id"].values})

        pm=df.groupby(["pid","cond"],observed=True)["x"].mean().reset_index()
        t=pm[pm.cond=="Tumor"]["x"]; b=pm[pm.cond=="Benign"]["x"]
        delta=float(t.mean()-b.mean()) if len(t) and len(b) else np.nan
        try: _,pv=mannwhitneyu(t,b) if len(t) and len(b) else (np.nan,np.nan)
        except Exception: pv=np.nan
        gene_rows.append(dict(cancer=tag,gene=g,delta=delta,tumor_mean=float(t.mean()) if len(t) else np.nan,
                              benign_mean=float(b.mean()) if len(b) else np.nan,p=pv))
gs=pd.DataFrame(gene_rows); gs.to_csv(os.path.join(TAB,"core_gene_tumor_induction.csv"),index=False)
print("=== per-gene tumor-induction (Δ tumor-benign, patient-level log-norm) ===")
print(gs.pivot_table(index="gene",columns="cancer",values="delta").reindex(CORE).to_string())
print("\nwrote figs: fig1_umap_core_score, fig2_core_detection, fig3_patient_boxplots, fig4_forest (300dpi)")
print("wrote core_gene_tumor_induction.csv")
