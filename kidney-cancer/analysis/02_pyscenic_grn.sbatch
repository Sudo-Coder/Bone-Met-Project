#!/bin/bash
#SBATCH --job-name=pyscenic
#SBATCH --account=ahussain-lab
#SBATCH --partition=all
#SBATCH --exclude=bellatrix,igsdgx-1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=10:00:00
#SBATCH --output=/autofs/projects-t3/hussain/scProj/kidney-cancer/outputs/logs/pyscenic-%j.out
#SBATCH --error=/autofs/projects-t3/hussain/scProj/kidney-cancer/outputs/logs/pyscenic-%j.err
set -euo pipefail

ROOT=/autofs/projects-t3/hussain/scProj
PROJ=$ROOT/kidney-cancer
PY=$ROOT/envs/scenic_env/bin
SCE=$PROJ/outputs/scenic
CT=$ROOT/resources/cistarget
N=${SLURM_CPUS_PER_TASK:-32}

# Prevent BLAS oversubscription — arboreto/dask spawns N single-thread workers.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1

LOOM=$SCE/myeloid_counts.loom
TFS=$CT/allTFs_hg38.txt
DB1=$CT/hg38_10kbp.genes_vs_motifs.rankings.feather
DB2=$CT/hg38_500bp.genes_vs_motifs.rankings.feather
ANNOT=$CT/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

echo "[$(date)] host=$(hostname) workers=$N"
cd $SCE

echo "[$(date)] STEP 1/3 GRNBoost2 ..."
$PY/pyscenic grn "$LOOM" "$TFS" -o adjacencies.tsv --num_workers $N --method grnboost2 --seed 777

echo "[$(date)] STEP 2/3 cisTarget prune ..."
$PY/pyscenic ctx adjacencies.tsv "$DB1" "$DB2" \
  --annotations_fname "$ANNOT" \
  --expression_mtx_fname "$LOOM" \
  --output regulons.csv \
  --num_workers $N --mask_dropouts

echo "[$(date)] STEP 3/3 AUCell ..."
$PY/pyscenic aucell "$LOOM" regulons.csv -o auc_mtx.csv --num_workers $N

echo "[$(date)] PYSCENIC DONE"
ls -lh adjacencies.tsv regulons.csv auc_mtx.csv
