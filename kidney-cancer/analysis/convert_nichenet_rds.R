#!/usr/bin/env Rscript
# Convert NicheNet v2 prior RDS files -> gzipped TSV so downstream NicheNet
# scoring can run in Python (no nichenetr / no modern R needed).
# Uses only base R (system R 3.6.3). Run:
#   Rscript analysis/convert_nichenet_rds.R
NN <- "/autofs/projects-t3/hussain/scProj/resources/nichenet"

cat("reading ligand_target_matrix...\n")
lt <- readRDS(file.path(NN, "ligand_target_matrix_nsga2r_final.rds"))
cat("  dim:", paste(dim(lt), collapse=" x "),
    "| rows(genes):", paste(head(rownames(lt),2), collapse=","),
    "| cols(ligands):", paste(head(colnames(lt),2), collapse=","), "\n")
# write as TSV with a leading 'target' column (rownames = target genes, cols = ligands)
gz <- gzfile(file.path(NN, "ligand_target_matrix.tsv.gz"), "w")
write.table(data.frame(target=rownames(lt), lt, check.names=FALSE),
            gz, sep="\t", quote=FALSE, row.names=FALSE)
close(gz)
cat("  wrote ligand_target_matrix.tsv.gz\n")

cat("reading lr_network...\n")
lr <- readRDS(file.path(NN, "lr_network_human_21122021.rds"))
cat("  class:", class(lr)[1], "| cols:", paste(colnames(lr), collapse=","), "| nrow:", nrow(lr), "\n")
write.table(as.data.frame(lr), gzfile(file.path(NN, "lr_network.tsv.gz"), "w"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("  wrote lr_network.tsv.gz\n")

cat("reading weighted_networks...\n")
wn <- readRDS(file.path(NN, "weighted_networks_nsga2r_final.rds"))
cat("  class:", class(wn)[1], "| names:", paste(names(wn), collapse=","), "\n")
if (is.list(wn)) {
  for (nm in names(wn)) {
    part <- as.data.frame(wn[[nm]])
    fn <- file.path(NN, paste0("weighted_networks_", nm, ".tsv.gz"))
    write.table(part, gzfile(fn, "w"), sep="\t", quote=FALSE, row.names=FALSE)
    cat("  wrote", basename(fn), "cols:", paste(colnames(part), collapse=","), "nrow:", nrow(part), "\n")
  }
}
cat("NICHENET CONVERT DONE\n")
