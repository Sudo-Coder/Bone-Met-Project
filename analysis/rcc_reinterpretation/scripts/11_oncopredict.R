#!/usr/bin/env Rscript
suppressMessages(library(oncoPredict))
root <- getwd()
gd <- file.path(root, "resources/model_data/oncopredict/gdsc1")
expr_train <- readRDS(file.path(gd, "GDSC1_Expr (RMA Normalized and Log Transformed).rds"))
res_train <- readRDS(file.path(gd, "GDSC1_Res.rds"))
res_train <- exp(res_train)
keep <- c("Sunitinib_5","Pazopanib_199","Axitinib_1021","Cabozantinib_249","Temsirolimus_1016")
res_train <- res_train[, keep]
kirc <- read.delim(gzfile(file.path(root, "resources/tcga/KIRC_HiSeqV2.gz")), row.names = 1, check.names = FALSE)
rownames(kirc) <- toupper(rownames(kirc))
tum <- colnames(kirc)[substr(colnames(kirc), 14, 15) %in% c("01","05")]
test_expr <- as.matrix(kirc[, tum])
od <- file.path(root, "analysis/rcc_reinterpretation/outputs/model/drug_out")
dir.create(od, showWarnings = FALSE, recursive = TRUE)
setwd(od)
pred <- calcPhenotype(trainingExprData = expr_train, trainingPtype = res_train, testExprData = test_expr,
              batchCorrect = "eb", powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2,
              minNumSamples = 10, printOutput = FALSE, removeLowVaringGenesFrom = "rawData")
cat("class:", class(pred), "dim:", paste(dim(pred), collapse = "x"), "\n")
if(!is.null(pred)) write.csv(pred, file.path(od, "DrugPredictions.csv"))
print(list.files(od, recursive = TRUE))
