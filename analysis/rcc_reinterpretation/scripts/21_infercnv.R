suppressMessages({library(infercnv); library(Matrix)})
set.seed(0)
ROOT <- "/autofs/projects-t3/hussain/scProj"
IN   <- file.path(ROOT, "analysis/rcc_reinterpretation/outputs/tables/copykat_input")
OUT  <- file.path(ROOT, "analysis/rcc_reinterpretation/outputs/cnv/infercnv")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cts   <- readMM(file.path(IN, "counts.mtx"))          # genes x cells
genes <- read.table(file.path(IN, "genes.txt"), sep = "\t")[, 1]
cells <- read.table(file.path(IN, "cells.txt"), sep = "\t")[, 1]
cts <- as(as(cts, "CsparseMatrix"), "dgCMatrix")   # readMM gives dgTMatrix; infercnv needs dgCMatrix
rownames(cts) <- genes; colnames(cts) <- cells
meta  <- read.csv(file.path(IN, "meta.csv"), stringsAsFactors = FALSE)
stopifnot(identical(meta$cell, cells))

ann <- data.frame(row.names = meta$cell, group = meta$compartment)
gof <- file.path(ROOT, "resources/annotation/gene_order_GRCh38.txt")

# T/NK are the reference. Myeloid is deliberately NOT a reference: whether TAMs are
# diploid is the thing being tested, so they must be scored as observations.
obj <- CreateInfercnvObject(raw_counts_matrix = cts,
                           annotations_file  = ann,
                           gene_order_file   = gof,
                           ref_group_names   = c("T_NK"))

obj <- infercnv::run(obj,
                     cutoff            = 0.1,      # 10x
                     out_dir           = OUT,
                     cluster_by_groups = TRUE,
                     denoise           = TRUE,
                     HMM               = TRUE,
                     HMM_type          = "i6",
                     analysis_mode     = "samples",
                     num_threads       = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")),
                     no_prelim_plot    = TRUE,
                     write_expr_matrix = TRUE,
                     save_rds          = TRUE)

# per-cell CNV burden = mean squared deviation of the denoised residuals from 1
resid <- read.table(file.path(OUT, "infercnv.observations.txt"), header = TRUE, check.names = FALSE)
ref   <- read.table(file.path(OUT, "infercnv.references.txt"),  header = TRUE, check.names = FALSE)
score <- function(m) colMeans((as.matrix(m) - 1)^2, na.rm = TRUE)
sc <- c(score(resid), score(ref))
df <- data.frame(cell = names(sc), cnv_score = as.numeric(sc))
df$cell <- gsub("\\.", "-", df$cell)
df <- merge(df, meta, by = "cell", all.x = TRUE)
write.csv(df, file.path(OUT, "cnv_score_per_cell.csv"), row.names = FALSE)

agg <- aggregate(cnv_score ~ compartment, df, function(x) c(n = length(x), mean = mean(x), median = median(x)))
print(do.call(data.frame, agg))
cat("[done] inferCNV\n")
