#!/bin/bash
#SBATCH --job-name=myeloid-build
#SBATCH --account=ahussain-lab
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:40:00
# CPU-only scanpy pipeline (HVG/scale/PCA/neighbors/leiden/umap); avoid slow/GPU nodes.
#SBATCH --exclude=bellatrix,igsdgx-1
#SBATCH --output=outputs/basic_analysis/logs/%x-%J.log
#SBATCH --error=outputs/basic_analysis/logs/%x-%J.err

ROOT=/autofs/projects-t3/hussain/scProj/kidney-cancer
ENV=/autofs/projects-t3/hussain/scProj/envs/scvi_gpu
cd "$ROOT" || exit 1
export PATH="$ENV/bin:$PATH"
export LD_LIBRARY_PATH="$ENV/lib:$LD_LIBRARY_PATH"
# keep BLAS within the allocation (scanpy steps are single-process vectorized)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "node: $(hostname -s)  start: $(date)"
"$ENV/bin/python" build_myeloid_final.py
echo "exit: $?  end: $(date)"
