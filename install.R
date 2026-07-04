# install.R — R environment for the RCC CLEC_LAM reinterpretation (R side of Phases 1–5).
# Mirrors the model papers' tooling: CellChat, (Multi)NicheNet, inferCNV, CopyKAT, GSVA/ssGSEA,
# UCell/AUCell, survival + survminer, ComplexHeatmap, plus AnnData<->Seurat/SCE interop.
# Run once on the login node (compute nodes have no internet). Prefer renv for a lockfile:
#   renv::init(); source("install.R"); renv::snapshot()  -> writes renv.lock
# Pin R >= 4.3. Set a scratch lib/cache if HOME is tight.

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

cran_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "data.table", "tidyverse",
  "survival", "survminer", "broom", "lme4", "lmerTest", "car",   # models + collinearity (VIF)
  "ggpubr", "patchwork", "pheatmap", "reticulate"
)
bioc_pkgs <- c(
  "GSVA",            # ssGSEA of CORE signature (mirrors Mei)
  "UCell", "AUCell", # rank-based single-cell scoring
  "infercnv",        # malignant sender calls (Phase 5)
  "ComplexHeatmap",
  "zellkonverter", "SingleCellExperiment", "anndata2ri", # .h5ad <-> SCE/Seurat
  "limma", "edgeR"   # pseudobulk DE if needed
)
github_pkgs <- c(
  "jinworks/CellChat",          # niche->TAM signaling (Phase 2)
  "saeyslab/nichenetr",         # ligand->target support (Phase 2)
  "saeyslab/multinichenetr",    # differential NicheNet (tumor vs benign)
  "navinlabcode/copykat"        # CNV cross-check (Phase 5)
)

install.packages(cran_pkgs)
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)
for (g in github_pkgs) try(remotes::install_github(g, upgrade = "never"))

# SeuratDisk / sceasy for object conversion (optional, alternative to zellkonverter)
try(remotes::install_github("mojaveazure/seurat-disk", upgrade = "never"))
try(remotes::install_github("cellgeni/sceasy", upgrade = "never"))

# NicheNet reference matrices are pre-staged at resources/nichenet/*.rds (no download on compute nodes).
# set.seed(0) in every script.
sessionInfo()
