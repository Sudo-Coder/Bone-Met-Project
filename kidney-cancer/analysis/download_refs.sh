#!/bin/bash
set -e
cd /autofs/projects-t3/hussain/scProj
CT=resources/cistarget; NN=resources/nichenet
BASE=https://resources.aertslab.org/cistarget
echo "[cistarget] rankings 10kb..."
curl -L --retry 3 -C - -o $CT/hg38_10kbp.genes_vs_motifs.rankings.feather \
  "$BASE/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
echo "[cistarget] rankings 500bp..."
curl -L --retry 3 -C - -o $CT/hg38_500bp.genes_vs_motifs.rankings.feather \
  "$BASE/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
echo "[cistarget] motif2tf + TF list..."
curl -L --retry 3 -o $CT/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl "$BASE/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
curl -L --retry 3 -o $CT/allTFs_hg38.txt "$BASE/tf_lists/allTFs_hg38.txt"
echo "[nichenet] priors..."
for f in ligand_target_matrix_nsga2r_final.rds lr_network_human_21122021.rds weighted_networks_nsga2r_final.rds; do
  curl -L --retry 3 -C - -o $NN/$f "https://zenodo.org/records/7074291/files/$f?download=1"
done
echo "DONE downloads"
ls -lh $CT $NN
