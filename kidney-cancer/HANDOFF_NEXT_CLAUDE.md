# Handoff prompt — paste this to the next Claude (scProj / kidney-cancer, mileena HPC)

You are picking up an scRNA-seq project on the IGS HPC (login node `mileena`, user `arifai`,
project at `/autofs/projects-t3/hussain/scProj`). Read your auto-memory (`MEMORY.md` and the
linked notes: `scproj-envs`, `scproj-layout-relocation`, `scvi-gpu-env`, `prostate-cancer-hpc`)
for full background. Here is the live state and the exact task.

## Background (why we're doing this)
The kidney-cancer myeloid CellOracle analysis needs to knock out **ATF3**, but ATF3 was missing
from the kidney data. Root cause (verified): the 9 **Benign-stroma** samples (`GSM6507*-Benign-stroma`)
were quantified against a reduced ~18k-gene annotation that LACKS ATF3 and the whole AP-1/immediate-early
module (JUN, JUNB, JUND, FOS, FOSB, EGR1, NR4A1, DUSP1). The old `basic_analysis` merged all samples
with `sc.concat(out)` (default `join='inner'` = gene INTERSECTION), so ATF3 got silently dropped.
Fix: drop the 9 stroma samples (they contribute 0 myeloid cells anyway) and use `join='outer'`.
The 17 RCC + 7 Benign-immune samples all carry the full 32,739-gene annotation WITH ATF3.

## What is running right now
A clean headless rerun of the kidney basic_analysis pipeline is on a V100:
- Script: `/autofs/projects-t3/hussain/scProj/kidney-cancer/pipeline_basic_analysis.py`
- Submitted via: `kidney-cancer/run_basic_analysis.sbatch`  (SLURM job **19154411**, node igsdgx-1)
- Check status: `export PATH=/usr/local/packages/slurm/bin:$PATH; squeue -j 19154411`
  Log: `kidney-cancer/outputs/basic_analysis/logs/basic-19154411.out` (look for `ALL DONE`).
- When done it writes (all in `kidney-cancer/Cleaned_Data/` unless noted):
  - `combined.h5ad`  (24 samples, outer join, ~32,739 genes, ATF3 present)
  - `integrated.h5ad`  (scVI-integrated, FULL gene set ~16k, ATF3 present; obs: `Sample`,
    `cell type` [from cell-annotations.csv, by barcode], `condition` [Tumor/Involved/Distal/Benign],
    `leiden`; obsm `X_scVI`; layers `counts` + `scvi_normalized`)
  - `kidney-cancer/scvi_integration_model/`  (saved scVI model)
  - `kidney-cancer/outputs/basic_analysis/{figures/*.png, leiden_markers.csv}`
- Originals preserved as `combined-with-stromal.h5ad` / `integrated-with-stromal.h5ad` (do NOT overwrite).

## YOUR TASK (do this AFTER job 19154411 finishes cleanly — confirm `ALL DONE` and `ATF3 in var_names`)
Build a NEW `kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad` using the SAME method the
**prostate** myeloid object was made, with a populated `final_label` column. The user will run
CellOracle themselves against this file.

Method (mirror prostate `single_cell_cellrank_myeloid_analysis.ipynb` c38-39 / celloracle c16):
1. Load the NEW `Cleaned_Data/integrated.h5ad`.
2. Subset to the myeloid cells = the barcodes present in `Cleaned_Data/myeloid-annotations.csv`
   (that CSV, `barcode,celltype`, already defines the 16,881 established myeloid cells and their
   labels; it was exported from the OLD myeloid object's `final_label`). Use `.reindex`/intersection
   on `obs_names` (barcodes are `RCC-BM10-Tumor_<bc>-1`, matching format — verified).
3. `myeloid.raw = myeloid`  (keep full gene set in .raw)
4. HVG subset like prostate: `sc.pp.highly_variable_genes(myeloid, subset=True, inplace=True)`
   (default dispersion flavor; prostate landed ~2,630 genes). **CRITICAL — before subsetting, FORCE
   ATF3 and the AP-1/IEG module into the kept set** so the whole point isn't lost to the top-HVG cut:
   set `highly_variable=True` for any of
   `['ATF3','JUN','JUNB','JUND','FOS','FOSB','EGR1','NR4A1','NR4A2','DUSP1']` present in `var_names`,
   THEN subset. Verify `'ATF3' in myeloid.var_names` after.
5. Recompute embedding on HVGs (prostate did): `sc.pp.scale(max_value=10)` → `sc.pp.pca` →
   `sc.pp.neighbors(use_rep='X_scVI' if you prefer the integrated latent, else PCA)` → `sc.tl.leiden`
   → `sc.tl.umap`. (Prostate used PCA-on-HVGs at n_neighbors=50,n_pcs=50,leiden res=0.5 — see
   cellrank c39. Match that unless the user says otherwise. Keep `X_scVI`, `counts`, `scvi_normalized`.)
6. Set `myeloid.obs['final_label']` from `myeloid-annotations.csv` by barcode:
   `final = pd.read_csv('Cleaned_Data/myeloid-annotations.csv', index_col='barcode')['celltype']`
   `myeloid.obs['final_label'] = pd.Categorical(final.reindex(myeloid.obs_names).fillna('Unknown').values)`
   (Labels: Mono1, Mono2, Mono3, TAM, TIM, mDC, Osteoclasts, Monocyte Pro. Report how many cells
   fell to 'Unknown' — the rerun's SOLO doublet calls differ slightly from the old object, so a few
   new cells won't match; that's expected.)
7. Keep `condition` (Tumor/Involved/Distal/Benign) and `Sample`. Ensure `counts` layer is present
   (CellOracle's `import_anndata_as_raw_count` uses raw counts + a `leiden` cluster column).
8. Write `Cleaned_Data/myeloid_FINAL_labels.h5ad`. Print final shape, `ATF3 present: True`, and the
   `final_label` value_counts.

Run it with the GPU env python (`/autofs/projects-t3/hussain/scProj/envs/scvi_gpu/bin/python`,
`torch.cuda` verified on V100). scVI isn't needed for this step (it's subset+HVG+cluster), so it can
run on the login node OR a short sbatch — your call; it's light. Prefer writing a small standalone
`.py` (do NOT edit the open notebooks — VSCode has been auto-reverting on-disk notebook edits;
scripts are safe).

## Gotchas / environment
- SLURM binaries not on PATH: `export PATH=/usr/local/packages/slurm/bin:$PATH`. Account = `ahussain-lab`
  (NEVER `igs`). GPU: `--gres=gpu:V100:1` lands on igsdgx-1 (the only V100 node).
- Big pip/torch installs must use scratch: `PIP_CACHE_DIR=/autofs/scratch/arifai/pip_cache`,
  `TMPDIR=/autofs/scratch/arifai/tmp` (HOME is only 10 GB).
- Compute nodes have NO internet (ribo list is pre-staged at `Cleaned_Data/ribo_genes_KEGG.txt`).
- Auto-approve is ON for this project (`.claude/settings.local.json` bypassPermissions + an allow-all
  PreToolUse hook that still hard-blocks catastrophic `rm -rf`/`sudo`/`mkfs`).
- The prostate project (`/autofs/projects-t3/hussain/scProj/prostate-cancer`) is the reference for
  "how it should look" — its `myeloid_FINAL.h5ad` is 2,630 genes with ATF3 (highly_variable=True).
