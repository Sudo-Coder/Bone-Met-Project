#!/bin/bash
#SBATCH --job-name=sc_ptime
#SBATCH --account=ahussain-lab
#SBATCH --partition=all
#SBATCH --exclude=bellatrix,igsdgx-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=/autofs/projects-t3/hussain/scProj/kidney-cancer/outputs/logs/sc_ptime-%j.out
#SBATCH --error=/autofs/projects-t3/hussain/scProj/kidney-cancer/outputs/logs/sc_ptime-%j.err
set -euo pipefail
ROOT=/autofs/projects-t3/hussain/scProj
PROJ=$ROOT/kidney-cancer
cd $PROJ
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "[$(date)] host=$(hostname)"
$ROOT/envs/cellrank_env/bin/python analysis/07_pseudotime_trends.py
echo "[$(date)] PTIME DONE"
