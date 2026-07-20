#!/usr/bin/env Rscript
lib <- path.expand("~/R_libs/4.4"); dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=max(1, parallel::detectCores()-2))
ok <- function(p) requireNamespace(p, quietly=TRUE)
inst_cran <- function(pk) for(p in pk) if(!ok(p)) tryCatch(install.packages(p, lib=lib), error=function(e) message("CRAN FAIL ",p,": ",conditionMessage(e)))
inst_bioc <- function(pk) for(p in pk) if(!ok(p)) tryCatch(BiocManager::install(p, lib=lib, update=FALSE, ask=FALSE), error=function(e) message("BIOC FAIL ",p,": ",conditionMessage(e)))
inst_gh   <- function(pk) for(p in pk) tryCatch(remotes::install_github(p, lib=lib, upgrade="never"), error=function(e) message("GH FAIL ",p,": ",conditionMessage(e)))

message("== CRAN deps ==")
inst_cran(c("NMF","ggalluvial","survminer","reticulate","anndata","future","future.apply","sna","ggnetwork","svglite","systemfonts"))
message("== Bioc ==")
inst_bioc(c("BiocNeighbors","GSVA","UCell","infercnv","zellkonverter","SingleCellExperiment","limma","edgeR","ComplexHeatmap"))
message("== GitHub: presto, CellChat, nichenetr, copykat ==")
inst_gh(c("immunogenomics/presto","jinworks/CellChat","saeyslab/nichenetr","navinlabcode/copykat"))

message("\n== VERIFY ==")
for(p in c("Seurat","CellChat","nichenetr","infercnv","copykat","GSVA","UCell","survminer","presto","NMF","zellkonverter"))
  cat(sprintf("  %-14s %s\n", p, if(ok(p)) "OK" else "MISSING"))
