#!/usr/bin/env Rscript
# 14_nichenet.R — Phase 2: ligand->target support for the CLEC_LAM receiver program (nichenetr).
# CellChat gives LR probability, not target causality. NicheNet tests whether sender-expressed ligands
# explain the CLEC_LAM-high transcriptional program (geneset_oi = up in TAM_CLEC_LAM vs TAM_other).
# Uses pre-staged resources/nichenet/*.rds (no download). Reports ranks of pre-specified ligands.
suppressMessages({library(nichenetr); library(Matrix); library(dplyr)})
lib <- path.expand("~/R_libs/4.4"); .libPaths(c(lib,.libPaths()))
ROOT <- getwd()
TAB <- file.path(ROOT,"analysis/rcc_reinterpretation/outputs/tables")
RES <- file.path(ROOT,"resources/nichenet")
ltm <- readRDS(file.path(RES,"ligand_target_matrix_nsga2r_final.rds"))
lr  <- readRDS(file.path(RES,"lr_network_human_21122021.rds"))
wn  <- readRDS(file.path(RES,"weighted_networks_nsga2r_final.rds"))
PRESPEC <- c("APOE","C1QA","C1QB","C1QC","C3","GAS6","PROS1","TGFB1")

d <- file.path(TAB,"cellchat_input/tumor")
expr <- as(readMM(file.path(d,"expr.mtx")),"CsparseMatrix")
rownames(expr) <- readLines(file.path(d,"genes.txt")); colnames(expr) <- readLines(file.path(d,"cells.txt"))
meta <- read.csv(file.path(d,"meta.csv"), row.names=1)
grp <- meta$cc_group; names(grp) <- rownames(meta)

frac_expr <- function(cells) { m <- expr[,cells,drop=FALSE]; Matrix::rowMeans(m>0) }
receiver <- names(grp)[grp=="TAM_CLEC_LAM"]
tam_other <- names(grp)[grp=="TAM_other"]
senders  <- names(grp)[grp %in% c("Tumor","MSC","Pericyte","Endothelial","Osteoclast","TAM_other","TAM_CLEC_LAM")]
expr_receiver <- frac_expr(receiver); expr_sender <- frac_expr(senders)
bg <- intersect(names(expr_receiver)[expr_receiver>0.10], rownames(ltm))
# geneset_oi = up in CLEC_LAM vs other TAM (mean log-norm diff), expressed, in ltm rownames
mu_c <- Matrix::rowMeans(expr[,receiver,drop=FALSE]); mu_o <- Matrix::rowMeans(expr[,tam_other,drop=FALSE])
diff <- mu_c - mu_o
geneset <- intersect(names(sort(diff,decreasing=TRUE))[1:300], bg)
geneset <- geneset[diff[geneset] > 0.25]
message("geneset_oi: ",length(geneset)," genes; background: ",length(bg))

# potential ligands = LR ligands expressed in senders whose receptor expressed in receiver
lig_all <- unique(lr$from); rec_all <- unique(lr$to)
expressed_lig <- intersect(lig_all, names(expr_sender)[expr_sender>0.10])
expressed_rec <- intersect(rec_all, names(expr_receiver)[expr_receiver>0.10])
lr_e <- lr %>% filter(from %in% expressed_lig & to %in% expressed_rec)
potential_ligands <- unique(lr_e$from); potential_ligands <- intersect(potential_ligands, colnames(ltm))
message("potential ligands: ",length(potential_ligands))

la <- predict_ligand_activities(geneset=geneset, background_expressed_genes=bg,
        ligand_target_matrix=ltm, potential_ligands=potential_ligands)
la <- la %>% arrange(desc(aupr_corrected)) %>% mutate(rank=row_number())
write.csv(la, file.path(TAB,"nichenet_ligand_activities.csv"), row.names=FALSE)

# ranks of pre-specified ligands
ps <- la %>% filter(test_ligand %in% PRESPEC) %>% select(test_ligand,aupr_corrected,pearson,rank)
ps$n_ligands <- nrow(la)
write.csv(ps, file.path(TAB,"nichenet_prespecified_ranks.csv"), row.names=FALSE)
cat("\n=== NicheNet pre-specified ligand ranks (of ",nrow(la)," ligands) ===\n"); print(as.data.frame(ps), row.names=FALSE)

# ligand-target links for pre-specified + top10
top_lig <- unique(c(intersect(PRESPEC,la$test_ligand), head(la$test_ligand,10)))
lt <- lapply(top_lig, function(L){
  tg <- intersect(geneset, colnames(ltm)); if(length(tg)==0) return(NULL)
  w <- ltm[L, tg]; data.frame(ligand=L, target=names(w), weight=as.numeric(w))
})
lt <- do.call(rbind, lt); lt <- lt[order(-lt$weight),]
write.csv(lt, file.path(TAB,"nichenet_ligand_target_links.csv"), row.names=FALSE)
cat("DONE nichenet\n")
