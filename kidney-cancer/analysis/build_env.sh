#!/bin/bash
set -e
cd /autofs/projects-t3/hussain/scProj
export MAMBA_ROOT_PREFIX=/autofs/projects-t3/hussain/scProj/.mamba
MM=./bin/micromamba
echo "[1/2] creating env python=3.10 ..."
$MM create -y -p envs/mechanism_env -c conda-forge python=3.10 pip "numpy<2" "pandas<2"
echo "[2/2] pip installing analysis stack ..."
envs/mechanism_env/bin/pip install --no-input \
  "pyscenic==0.12.1" "decoupler" "omnipath" "scanpy" "pyarrow" "gseapy" \
  "statsmodels" "adjustText" "leidenalg" "python-igraph" "lifelines"
echo "DONE build_env"
envs/mechanism_env/bin/python -c "import pyscenic, decoupler, scanpy, omnipath; print('pyscenic', pyscenic.__version__); print('decoupler', decoupler.__version__); print('scanpy', scanpy.__version__)"
