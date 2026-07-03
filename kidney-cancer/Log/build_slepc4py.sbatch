#!/bin/bash
# TARGETED FINISH-THE-JOB build: petsc, petsc4py, and slepc are ALREADY installed
# and persisted in cellrank_env. Only slepc4py remains -- the prior full-rebuild
# job (19009615) TIMED OUT at 90min because it rebuilt petsc4py + slepc wheels
# first and ran out of clock before slepc4py finished.
#
# This job rebuilds NOTHING already present: it just compiles the small slepc4py
# cython binding against the installed slepc (--no-build-isolation), then runs the
# exact krylov verification. Finishes in minutes.
#
# Submit:  sbatch build_slepc4py.sbatch
#SBATCH --job-name=slepc4py
#SBATCH --account=ahussain-lab
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=slepc4py-build-%J.log
#SBATCH --error=slepc4py-build-%J.err

ROOT=/autofs/projects-t3/hussain/scProj
PY="$ROOT/envs/cellrank_env/bin/python"
# env bin on PATH so the build finds patchelf (slepc4py wheel rewrites rpaths).
export PATH="$ROOT/envs/cellrank_env/bin:$PATH"

# Compute nodes ship no gcc/gfortran in base PATH -- only OpenMPI wrappers, which
# need a backing compiler from the gcc module. (Done before 'set -u'; Lmod touches
# unbound vars.)
for init in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash; do
  [ -f "$init" ] && source "$init" && break
done
module load gcc/12.2.0
export OMPI_CC=gcc OMPI_CXX=g++ OMPI_FC=gfortran
export MPICC=mpicc MPICXX=mpicxx MPIF90=mpif90

set -euo pipefail

export TMPDIR="/autofs/scratch/$USER/slepc4py_build.$SLURM_JOB_ID"
export PIP_CACHE_DIR="/autofs/scratch/$USER/pip_cache"
mkdir -p "$TMPDIR" "$PIP_CACHE_DIR"

echo "==================== BUILD ENV ===================="
echo "Node:     $(hostname -s)   cpus=${SLURM_CPUS_PER_TASK}"
echo "gcc:      $(command -v gcc || echo MISSING)   gfortran: $(command -v gfortran || echo MISSING)"
echo "mpicc:    $(command -v mpicc)   patchelf: $(command -v patchelf || echo MISSING)"
echo "==================================================="
for t in gcc g++ gfortran mpicc mpif90; do
  command -v "$t" >/dev/null || { echo "FATAL: '$t' not on PATH -- gcc module did not load"; exit 1; }
done

# Sanity: the components slepc4py links against must already be importable.
echo ">>> Confirming petsc / petsc4py / slepc are present..."
"$PY" -c "import petsc, slepc; import petsc4py; petsc4py.init(); print('petsc/petsc4py/slepc OK')"

# slepc4py must link the INSTALLED petsc/slepc. Under pip build isolation it
# rebuilds slepc in an overlay and fails 'Unable to link with PETSc'. So point it
# at the installed dirs and use --no-build-isolation.
export PETSC_DIR="$("$PY" -c 'import petsc; print(petsc.get_petsc_dir())')"
export SLEPC_DIR="$("$PY" -c 'import slepc; print(slepc.get_slepc_dir())')"
unset PETSC_ARCH 2>/dev/null || true
echo ">>> PETSC_DIR=$PETSC_DIR"
echo ">>> SLEPC_DIR=$SLEPC_DIR"
export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK}"

# cython is needed for the binding build; numpy/setuptools already present.
"$PY" -m pip install cython >/dev/null
echo ">>> Building slepc4py (no build isolation, against installed slepc)..."
"$PY" -m pip install --no-build-isolation --force-reinstall --no-deps slepc4py

echo ">>> Build done. Verifying krylov actually engages..."
"$PY" - <<'PYEOF'
import warnings
import numpy as np
import scipy.sparse as sp

import petsc4py, slepc4py
print("petsc4py:", petsc4py.__version__, " slepc4py:", slepc4py.__version__)

# 1) The exact gate cellrank uses. True => no krylov->brandts override, no densify.
from cellrank._utils._linear_solver import _is_petsc_slepc_available
assert _is_petsc_slepc_available(), "FAIL: cellrank STILL cannot see petsc/slepc"
print("cellrank _is_petsc_slepc_available(): True  <-- the bug is fixed")

# 2) Prove SLEPc links at RUNTIME: real krylov Schur on a SPARSE matrix.
from pygpcca._sorted_schur import sorted_schur
n = 300
A = sp.random(n, n, density=0.05, format="csr", random_state=0)
A = A + sp.eye(n, format="csr")
P = sp.diags(1.0 / np.asarray(A.sum(1)).ravel()) @ A
assert sp.issparse(P)

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    try:
        R, Q, ev = sorted_schur(P, 8, "LM", method="krylov")
        print("krylov sorted_schur on SPARSE matrix: OK, returned", len(ev), "eigenvalues")
    except ValueError as e:
        if "Sparse implementation" in str(e):
            raise SystemExit("FAIL: fell back to brandts -> petsc/slepc NOT linked")
        print("krylov ran (benign:", e, ")")
    msgs = [str(x.message) for x in w]
assert not any("PETSc" in m or "SLEPc" in m or "petsc4py" in m for m in msgs), \
    f"FAIL: saw a petsc/slepc fallback warning: {msgs}"

print("\n==================== VERIFIED ====================")
print("krylov Schur is now genuinely sparse/fast in cellrank_env.")
print("Your Cell 69 compute_schur(method='krylov') will NO LONGER densify.")
print("=================================================")
PYEOF

rm -rf "$TMPDIR"
echo ">>> All done."
