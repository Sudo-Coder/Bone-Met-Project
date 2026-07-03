#!/usr/bin/env python
"""
Kidney-cancer basic_analysis pipeline — headless, GPU scVI, end-to-end up to
(but NOT including) the COVID/lung exploratory section of the notebook.

Faithful to single_cell_basic_analysis.ipynb's REAL pipeline cells, with the
fixes we agreed on:
  * 24 samples: 17 RCC + 7 Benign-immune, NO Benign-stroma  (stroma lacks ATF3)
  * sc.concat(..., join='outer')  -> ATF3 + AP-1 genes retained
  * labels applied by BARCODE from Cleaned_Data/cell-annotations.csv (the
    notebook's leiden->dict at c106 was empty), so labels survive new clustering
  * ribosomal gene list read from a LOCAL file (no internet on compute nodes)
  * matplotlib Agg + figures saved to outputs/ (viewable later in VSCode)

Stages (each writes a checkpoint; set REBUILD_COMBINED=1 to force stage 1):
  1. per-sample QC + SOLO doublet removal -> Cleaned_Data/combined.h5ad
  2. scVI integration + leiden + barcode labels -> Cleaned_Data/integrated.h5ad

Run via run_basic_analysis.sbatch on a V100.
"""
import os, sys, time
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import matplotlib
matplotlib.use("Agg")
import scanpy as sc
import scvi

# ----------------------------------------------------------------------------
PROJ = "/autofs/projects-t3/hussain/scProj/kidney-cancer"
os.chdir(PROJ)
CLEAN = f"{PROJ}/Cleaned_Data"
DATA  = f"{PROJ}/Data"
OUT   = f"{PROJ}/outputs/basic_analysis"
FIG   = f"{OUT}/figures"
os.makedirs(FIG, exist_ok=True)
sc.settings.figdir = FIG
sc.settings.autoshow = False
scvi.settings.seed = 0          # reproducibility (notebook was unseeded)

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

# 24 samples: RCC + Benign-immune only (stroma dropped)
FILES = [
 'GSM6133738_RCC-BM10-Tumor.count.csv', 'GSM6133739_RCC-BM10-Involve.count.csv',
 'GSM6133740_RCC-BM10-Noninvolved.count.csv', 'GSM6133741_RCC-BM3-Tumor.count.csv',
 'GSM6133742_RCC-BM4-Tumor.count.csv', 'GSM6133743_RCC-BM5-Tumor.count.csv',
 'GSM6133744_RCC-BM7-Tumor.count.csv', 'GSM6133745_RCC-BM8-Tumor.count.csv',
 'GSM6133746_RCC-BM1-Tumor.count.csv', 'GSM6133747_RCC-BM1-Involve.count.csv',
 'GSM6133748_RCC-BM1-Noninvolved.count.csv', 'GSM6133749_RCC-BM2-Tumor.count.csv',
 'GSM6133750_RCC-BM2-Involve.count.csv', 'GSM6133751_RCC-BM2-Noninvolved.count.csv',
 'GSM6133752_RCC-BM9-Involve.count.csv', 'GSM6133753_RCC-BM9-Noninvolved.count.csv',
 'GSM6133754_RCC-BM9-Tumor.count.csv', 'GSM6507012_BMM2-Benign-immune.count.csv',
 'GSM6507013_BMM3-Benign-immune.count.csv', 'GSM6507014_BMM4-Benign-immune.count.csv',
 'GSM6507015_BMM5-Benign-immune.count.csv', 'GSM6507016_BMM6-Benign-immune.count.csv',
 'GSM6507017_BMM8-Benign-immune.count.csv', 'GSM6507018_BMM9-Benign-immune.count.csv',
]

# ribosomal gene list, pre-downloaded (KEGG_RIBOSOME): line1=name, line2=url, rest=genes
ribo_genes = pd.read_csv(f"{CLEAN}/ribo_genes_KEGG.txt", skiprows=2, header=None)
ribo_gene_list = ribo_genes[0].values

# ---------------------------------------------------------------------------
def pp(csv_path):
    """Per-sample: HVG->scVI->SOLO doublet call, then QC-filter full-gene object.
    Faithful to notebook cell 70. scVI/SOLO train on the GPU (accelerator='auto')."""
    adata = sc.read_csv(csv_path).T
    adata.var_names_make_unique(); adata.obs_names_make_unique()
    adata.var['ribo'] = adata.var_names.isin(ribo_gene_list)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor='seurat_v3')
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata); vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae); solo.train()
    df = solo.predict(); df['prediction'] = solo.predict(soft=False)
    df.index = df.index.map(lambda x: x[:-2])
    df['dif'] = df.doublet - df.singlet
    doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]

    adata = sc.read_csv(csv_path).T
    adata.obs['Sample'] = os.path.basename(csv_path).replace('.csv', '')
    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]
    adata.obs_names_make_unique(); adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    # NOTE: original used lowercase 'mt-' (human mito is 'MT-'), so pct_counts_mt
    # is ~0 and the mt<20 filter is effectively inactive. Kept AS-IS for fidelity
    # to the object you're reproducing. Change to 'MT-' if you want real mt QC.
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribo'] = adata.var_names.isin(ribo_gene_list)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None,
                               log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata = adata[adata.obs.pct_counts_mt < 20]
    adata = adata[adata.obs.pct_counts_ribo < 100]
    adata.X = csr_matrix(adata.X)
    return adata

# ---------------------------------------------------------------------------
def stage1_combined():
    ckpt = f"{CLEAN}/combined.h5ad"
    if os.path.exists(ckpt) and os.environ.get("REBUILD_COMBINED") != "1":
        log(f"stage 1 SKIP (exists): {ckpt}  (set REBUILD_COMBINED=1 to force)")
        return
    out = []
    for i, f in enumerate(FILES, 1):
        log(f"stage 1 [{i}/{len(FILES)}] pp({f})")
        out.append(pp(os.path.join(DATA, f)))
    log("stage 1 concat (join='outer')")
    adata = sc.concat(out, join='outer')
    sc.pp.filter_genes(adata, min_cells=3)
    adata.X = csr_matrix(adata.X)
    adata.write_h5ad(ckpt)
    log(f"stage 1 DONE -> {ckpt}  shape={adata.shape}  ATF3 in genes: {'ATF3' in adata.var_names}")

# ---------------------------------------------------------------------------
def stage2_integrate():
    ckpt = f"{CLEAN}/integrated.h5ad"
    log("stage 2 read combined")
    adata = sc.read_h5ad(f"{CLEAN}/combined.h5ad")
    sc.pp.filter_genes(adata, min_cells=100)
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    log("stage 2 scVI integration (GPU)")
    scvi.model.SCVI.setup_anndata(
        adata, layer="counts", categorical_covariate_keys=["Sample"],
        continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'])
    model = scvi.model.SCVI(adata)
    model.train()
    adata.obsm['X_scVI'] = model.get_latent_representation()
    adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4)

    log("stage 2 save scVI integration model")
    model.save(f"{PROJ}/scvi_integration_model", overwrite=True)

    log("stage 2 neighbors / umap / leiden(res=0.9)")
    sc.pp.neighbors(adata, use_rep='X_scVI')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.9)
    sc.tl.rank_genes_groups(adata, 'leiden')
    markers = sc.get.rank_genes_groups_df(adata, None)
    markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
    adata.uns['markers'] = markers
    markers.to_csv(f"{OUT}/leiden_markers.csv", index=False)

    # ---- LABELS by barcode from cell-annotations.csv (broad atlas labels) ----
    log("stage 2 apply labels by barcode from cell-annotations.csv")
    annot = pd.read_csv(f"{CLEAN}/cell-annotations.csv", index_col='barcode')
    lab = annot['celltype'].reindex(adata.obs_names)
    n_missing = int(lab.isna().sum())
    lab = lab.fillna('Unknown')
    adata.obs['cell type'] = pd.Categorical(lab.values)
    log(f"stage 2 labels: {n_missing}/{adata.n_obs} cells had no annotation -> 'Unknown'")

    # ---- corrected kidney condition from Sample (replaces COVID map_condition) ----
    def kidney_condition(s):
        if 'Tumor' in s:        return 'Tumor'
        if 'Involve' in s:      return 'Involved'
        if 'Noninvolved' in s:  return 'Distal'
        if 'Benign' in s:       return 'Benign'
        return 'NA'
    adata.obs['condition'] = pd.Categorical(adata.obs['Sample'].map(kidney_condition))

    # ---- figures (saved, not shown) ----
    for color, fn in [(['leiden'], 'umap_leiden.png'),
                      (['Sample'], 'umap_sample.png'),
                      (['cell type'], 'umap_celltype.png'),
                      (['condition'], 'umap_condition.png')]:
        try:
            sc.pl.umap(adata, color=color, frameon=False, show=False, save=f"_{fn}")
        except Exception as e:
            log(f"  (figure {fn} failed: {e})")

    adata.write_h5ad(ckpt)
    log(f"stage 2 DONE -> {ckpt}  shape={adata.shape}  ATF3 present: {'ATF3' in adata.var_names}")

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    t0 = time.time()
    log(f"START  cuda_available={__import__('torch').cuda.is_available()}")
    stage1_combined()
    stage2_integrate()
    log(f"ALL DONE in {(time.time()-t0)/60:.1f} min")
