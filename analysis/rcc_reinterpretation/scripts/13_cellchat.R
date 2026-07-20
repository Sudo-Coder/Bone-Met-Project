#!/usr/bin/env Rscript
suppressMessages({library(CellChat); library(Matrix)})
lib <- path.expand("~/R_libs/4.4"); .libPaths(c(lib,.libPaths()))
ROOT <- normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),"..","..",".."))
if(length(ROOT)==0 || is.na(ROOT)) ROOT <- getwd()
TAB <- file.path(ROOT,"analysis/rcc_reinterpretation/outputs/tables")
IND <- file.path(TAB,"cellchat_input")
RECV <- c("TAM_CLEC_LAM","TAM_other","TIM")
PRESPEC <- c("COMPLEMENT","GAS","PROS","TGFb","SPP1","MIF","GALECTIN","ANNEXIN","GRN","APP")
future::plan("sequential")
options(stringsAsFactors=FALSE)

run_niche <- function(niche){
  d <- file.path(IND,niche)
  if(!file.exists(file.path(d,"expr.mtx"))){ message("skip ",niche); return(invisible()) }
  expr <- as(readMM(file.path(d,"expr.mtx")),"CsparseMatrix")
  rownames(expr) <- readLines(file.path(d,"genes.txt"))
  colnames(expr) <- readLines(file.path(d,"cells.txt"))
  meta <- read.csv(file.path(d,"meta.csv"), row.names=1)
  meta$cc_group <- factor(meta$cc_group)
  keep <- names(which(table(meta$cc_group) >= 10))
  cells <- rownames(meta)[meta$cc_group %in% keep]
  expr <- expr[,cells]; meta <- meta[cells,,drop=FALSE]; meta$cc_group <- droplevels(meta$cc_group)
  message(niche," groups: ", paste(levels(meta$cc_group),collapse=", "))
  cc <- createCellChat(object=expr, meta=meta, group.by="cc_group")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type="triMean")
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  saveRDS(cc, file.path(TAB,paste0("cellchat_",niche,".rds")))
  df <- subsetCommunication(cc)                     # LR-level data.frame
  write.csv(df, file.path(TAB,paste0("cellchat_allLR_",niche,".csv")), row.names=FALSE)
  tam <- df[df$target %in% RECV, ]
  tam <- tam[order(-tam$prob), ]
  write.csv(tam, file.path(TAB,paste0("cellchat_tam_LR_",niche,".csv")), row.names=FALSE)
  ax <- tam[tam$pathway_name %in% PRESPEC, ]
  write.csv(ax, file.path(TAB,paste0("cellchat_prespecified_axes_",niche,".csv")), row.names=FALSE)
  ct <- cc@netP$centr
  cen <- do.call(rbind, lapply(names(ct), function(pw){
    z <- ct[[pw]]; data.frame(pathway=pw, group=names(z$outdeg), outdeg=z$outdeg, indeg=z$indeg)
  }))
  write.csv(cen, file.path(TAB,paste0("cellchat_centrality_",niche,".csv")), row.names=FALSE)
  message(niche,": ",nrow(df)," LR pairs, ",nrow(tam)," to TAM receivers, ",nrow(ax)," pre-specified.")
  cat("TOP sender->TAM_CLEC_LAM axes:\n")
  x <- tam[tam$target=="TAM_CLEC_LAM",][1:min(15,sum(tam$target=='TAM_CLEC_LAM')),
             c("source","target","ligand","receptor","pathway_name","prob","pval")]
  print(x, row.names=FALSE)
}
for(n in c("tumor","benign")) tryCatch(run_niche(n), error=function(e) message("ERR ",n,": ",conditionMessage(e)))
cat("DONE cellchat\n")
