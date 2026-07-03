#!/bin/bash
set -e
cd /autofs/projects-t3/hussain/scProj
export MAMBA_ROOT_PREFIX=/autofs/projects-t3/hussain/scProj/.mamba
MM=./bin/micromamba
echo "[1/2] create scenic_env python=3.10"
$MM create -y -p envs/scenic_env -c conda-forge python=3.10 pip
echo "[2/2] pip install pinned pyscenic stack (numpy<1.24, pandas<2)"
envs/scenic_env/bin/pip install --no-input \
  "setuptools<81" "numpy==1.23.5" "pandas==1.5.3" "scipy==1.10.1" \
  "numba==0.56.4" "pyarrow==11.0.0" "dask==2023.3.2" "distributed==2023.3.2" \
  "pyscenic==0.12.1" "loompy==3.0.8"
echo "SCENIC ENV DONE"
envs/scenic_env/bin/python -c "import pyscenic,arboreto,ctxcore,numpy,pandas; print('pyscenic',pyscenic.__version__,'numpy',numpy.__version__,'pandas',pandas.__version__)"
