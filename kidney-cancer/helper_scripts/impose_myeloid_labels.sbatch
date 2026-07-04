#!/bin/bash
#SBATCH --job-name=impose-labels
#SBATCH --account=ahussain-lab
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:25:00
#SBATCH --exclude=bellatrix,igsdgx-1
#SBATCH --output=outputs/basic_analysis/logs/%x-%J.log
#SBATCH --error=outputs/basic_analysis/logs/%x-%J.err
ROOT=/autofs/projects-t3/hussain/scProj/kidney-cancer
ENV=/autofs/projects-t3/hussain/scProj/envs/sc-general-analysis
cd "$ROOT" || exit 1
export PATH="$ENV/bin:$PATH"
echo "node: $(hostname -s) start: $(date)"
"$ENV/bin/python" impose_myeloid_labels.py
echo "exit: $? end: $(date)"
