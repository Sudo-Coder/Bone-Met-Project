#!/bin/bash
#SBATCH --job-name=sc_cpu
#SBATCH --account=ahussain-lab
#SBATCH --partition=all
#SBATCH --exclude=bellatrix,igsdgx-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=06:00:00
#SBATCH --output=/autofs/projects-t3/hussain/scProj/kidney-cancer/outputs/logs/sc_cpu-%j.out
#SBATCH --error=/autofs/projects-t3/hussain/scProj/kidney-cancer/outputs/logs/sc_cpu-%j.err
set -euo pipefail
ROOT=/autofs/projects-t3/hussain/scProj
PROJ=$ROOT/kidney-cancer
PY=$ROOT/envs/mechanism_env/bin/python
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd $PROJ
echo "[$(date)] host=$(hostname)"
echo "[$(date)] === Task 3: decoupleR ==="
$PY analysis/04_decoupler.py
echo "[$(date)] === Task 2: NicheNet ==="
$PY analysis/06_nichenet.py
echo "[$(date)] CPU TASKS DONE"
