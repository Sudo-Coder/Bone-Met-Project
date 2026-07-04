#!/usr/bin/env Rscript
# 20b_copykat.R — Phase 1.5: CopyKAT malignant-sender validation (RCC full-niche subset).
# Calls aneuploid (malignant) vs diploid per cell using T/NK as known-normal reference.
# Expect: Tumor-labeled cells aneuploid; MSC/Pericyte/Endothelial/Myeloid/T diploid -> justifies
# the CellChat malignant sender = Tumor, non-malignant stromal senders.
suppressMessages({library(copykat); library(Matrix)})
lib <- path.expand("~/R_libs/4.4"); .libPaths(c(lib,.libPaths()))
ROOT <- getwd(); TAB <- file.path(ROOT,"analysis/rcc_reinterpretation/outputs/tables")
d <- file.path(TAB,"copykat_input")
m <- readMM(file.path(d,"counts.mtx")); rownames(m)<-readLines(file.path(d,"genes.txt")); colnames(m)<-readLines(file.path(d,"cells.txt"))
norm <- readLines(file.path(d,"normal_cells.txt"))
meta <- read.csv(file.path(d,"meta.csv"))
setwd(TAB)  # copykat writes outputs to cwd
res <- copykat(rawmat=as.matrix(m), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1,
               sam.name="rcc_ck", distance="euclidean", norm.cell.names=norm,
               n.cores=max(1,parallel::detectCores()-2), output.seg=FALSE, plot.genes=FALSE)
pred <- as.data.frame(res$prediction)
mm <- merge(pred, meta, by.x="cell.names", by.y="cell", all.x=TRUE)
write.csv(mm, file.path(TAB,"copykat_predictions.csv"), row.names=FALSE)
tab <- table(mm$compartment, mm$copykat.pred)
write.csv(as.data.frame.matrix(tab), file.path(TAB,"copykat_compartment_x_pred.csv"))
cat("\n=== CopyKAT: compartment x prediction ===\n"); print(tab)
cat("\nTumor aneuploid fraction:", round(mean(mm$copykat.pred[mm$compartment=='Tumor']=='aneuploid',na.rm=TRUE),3),"\n")
cat("Stroma (MSC+Peri+Endo) aneuploid fraction:",
    round(mean(mm$copykat.pred[mm$compartment %in% c('MSC','Pericyte','Endothelial')]=='aneuploid',na.rm=TRUE),3),"\n")
cat("DONE copykat\n")
