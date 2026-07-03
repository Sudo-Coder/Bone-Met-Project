#!/bin/bash
# Headless kidney basic_analysis pipeline on a V100.
# Submit:   export PATH=/usr/local/packages/slurm/bin:$PATH
#           sbatch kidney-cancer/run_basic_analysis.sbatch
# Watch:    tail -f kidney-cancer/outputs/basic_analysis/logs/basic-<JOBID>.out
#SBATCH --job-name=kid-basic
#SBATCH --account=ahussain-lab
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:V100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --exclude=bellatrix
#SBATCH --output=/autofs/projects-t3/hussain/scProj/kidney-cancer/outputs/basic_analysis/logs/basic-%J.out
#SBATCH --error=/autofs/projects-t3/hussain/scProj/kidney-cancer/outputs/basic_analysis/logs/basic-%J.err
set -euo pipefail

ROOT=/autofs/projects-t3/hussain/scProj
PROJ=$ROOT/kidney-cancer
PY=$ROOT/envs/scvi_gpu/bin/python

# keep torch/tmp off the tiny HOME + /tmp
export TMPDIR=/autofs/scratch/$USER/tmp
mkdir -p "$TMPDIR"
# BLAS threading for the CPU-side scanpy steps (neighbors/umap/leiden); scVI is on GPU
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK \
       OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd "$PROJ"
echo "[$(date)] host=$(hostname) gpu=${CUDA_VISIBLE_DEVICES:-none}"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# To force a full rebuild of combined.h5ad (re-run all 24 SOLO doublet calls):
#   export REBUILD_COMBINED=1
"$PY" "$PROJ/pipeline_basic_analysis.py"
echo "[$(date)] pipeline exited $?"
