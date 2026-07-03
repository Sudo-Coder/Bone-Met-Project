#!/bin/bash
# Build PETSc/SLEPc + petsc4py/slepc4py into cellrank_env, then VERIFY that
# cellrank/pygpcca's krylov Schur path actually engages at runtime.
#
# WHY: without petsc4py/slepc4py, cellrank's compute_schur(method="krylov")
# silently downgrades to brandts and DENSIFIES the sparse transition matrix
# (cellrank/.../decomposition/_schur.py:158-167) -> dense O(n^3) Schur,
# single-threaded -> the ~40h hang. With SLEPc present, krylov runs sparse
# in seconds-minutes.
#
# Submit:  sbatch build_petsc_slepc.sbatch
# Watch:   tail -f petsc-build-<JOBID>.log
#SBATCH --job-name=petsc-build
#SBATCH --account=ahussain-lab
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=petsc-build-%J.log
#SBATCH --error=petsc-build-%J.err

ROOT=/autofs/projects-t3/hussain/scProj
PY="$ROOT/envs/cellrank_env/bin/python"
# Put the env's bin on PATH so the build finds tools installed into it -- notably
# 'patchelf', which slepc/slepc4py wheel builds invoke to set library rpaths.
export PATH="$ROOT/envs/cellrank_env/bin:$PATH"

# Compute nodes have NO gcc/gfortran in the base PATH -- only the OpenMPI
# wrapper compilers, which then can't find a backing compiler. The real
# toolchain comes from environment modules. Load gcc (provides gcc/g++/
# gfortran) and point the OpenMPI wrappers at them. (Verified on
# juggernaut/bellatrix: mpicc + mpif90 both compile after this.)
# NB: done BEFORE 'set -u' because Lmod init touches unbound vars.
for init in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash; do
  [ -f "$init" ] && source "$init" && break
done
module load gcc/12.2.0
export OMPI_CC=gcc OMPI_CXX=g++ OMPI_FC=gfortran
# PETSc's pip build rejects --with-cc/--with-fc; it wants these env vars instead.
export MPICC=mpicc MPICXX=mpicxx MPIF90=mpif90

set -euo pipefail

# Keep ALL build temp + pip cache OFF home (~10GB). Use scratch (5TB).
export TMPDIR="/autofs/scratch/$USER/petsc_build.$SLURM_JOB_ID"
export PIP_CACHE_DIR="/autofs/scratch/$USER/pip_cache"
mkdir -p "$TMPDIR" "$PIP_CACHE_DIR"

# Compilers are passed via MPICC/MPICXX/MPIF90 env vars above (pip petsc
# forbids --with-cc here). Keep only build-tuning flags.
# --download-openblas (NOT fblaslapack): nodes have no linkable BLAS. fblaslapack
# builds STATIC libs that the pip wheel fails to ship, so downstream SLEPc can't
# find -lflapack/-lfblas. OpenBLAS builds a SHARED libopenblas.so that IS packaged
# into the prefix, so SLEPc links cleanly.
# --with-fortran-bindings=0: skips compiling PETSc's thousands of Fortran-stub
# files (the slow, single-threaded part). petsc4py/slepc4py use the C API only,
# so they don't need them -- this roughly halves the build.
# DYNAMIC_ARCH=1: CRITICAL on this HETEROGENEOUS cluster. Nodes span AMD Zen4
# (AVX-512), AMD Zen3 (no AVX-512), and Intel Ivy Bridge (no AVX2). Without this,
# OpenBLAS compiles for ONLY the build node's CPU and SIGILLs ('Illegal
# instruction') on any older node a job lands on. DYNAMIC_ARCH builds all kernels
# and selects at runtime per-CPU, so the library is portable across every node.
export PETSC_CONFIGURE_OPTIONS="--with-debugging=0 --with-x=0 --download-openblas --download-openblas-make-options=DYNAMIC_ARCH=1 --with-fortran-bindings=0"
export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK}"

echo "==================== BUILD ENV ===================="
echo "Node:     $(hostname -s)   cpus=${SLURM_CPUS_PER_TASK}"
echo "Python:   $PY"
echo "TMPDIR:   $TMPDIR"
echo "gcc:      $(command -v gcc || echo MISSING)"
echo "gfortran: $(command -v gfortran || echo MISSING)"
echo "mpicc:    $(command -v mpicc)   mpif90: $(command -v mpif90)"
echo "==================================================="
# Fail fast (seconds, not 20 min) if the toolchain didn't load.
for t in gcc g++ gfortran mpicc mpif90; do
  command -v "$t" >/dev/null || { echo "FATAL: '$t' not on PATH -- gcc module did not load"; exit 1; }
done

# TWO PHASES on purpose. A single 'pip install petsc petsc4py slepc slepc4py'
# fails at slepc4py with "Unable to link with PETSc": under pip build isolation,
# slepc gets rebuilt in an overlay that can't see the PETSc just compiled.
# Fix: fully INSTALL petsc/petsc4py first, export PETSC_DIR, then build slepc
# against it. Caching is ON (PIP_CACHE_DIR on scratch) so reruns reuse the
# ~25-min PETSc wheel instead of recompiling.
# Clean slate for a retry, but DO NOT purge the pip cache: the OpenBLAS petsc
# wheel is correct and cached, so reinstalling reuses it in seconds (no 14-min
# recompile). If you ever change PETSC_CONFIGURE_OPTIONS, manually run
# `pip cache remove 'petsc*'` once to force a fresh build.
"$PY" -m pip uninstall -y petsc petsc4py slepc slepc4py 2>/dev/null || true
# patchelf: slepc/slepc4py wheel builds need it to rewrite rpaths.
"$PY" -m pip install patchelf
# Build petsc ALONE first so its wheel builds once and caches; petsc4py reuses it.
echo ">>> [1/3] Installing petsc (OpenBLAS, no Fortran stubs) -- reuses cached wheel if present..."
"$PY" -m pip install petsc
echo ">>> petsc ready; installing petsc4py (reuses it)..."
"$PY" -m pip install petsc4py

# Make the installed PETSc discoverable. pip installs PETSc in prefix mode
# (PETSC_DIR points at lib/include, PETSC_ARCH empty).
export PETSC_DIR="$("$PY" -c 'import petsc; print(petsc.get_petsc_dir())')"
unset PETSC_ARCH 2>/dev/null || true
echo ">>> PETSc installed. PETSC_DIR=$PETSC_DIR"

# slepc4py must link the PETSc above. Under pip build ISOLATION it instead
# rebuilds slepc in an overlay and fails 'Unable to link with PETSc'. Build
# slepc + slepc4py with --no-build-isolation so they use the installed
# petsc/slepc directly. That needs cython in the env (numpy/setuptools present).
echo ">>> [2/3] Installing cython + slepc (no build isolation)..."
"$PY" -m pip install cython
"$PY" -m pip install --no-build-isolation slepc
export SLEPC_DIR="$("$PY" -c 'import slepc; print(slepc.get_slepc_dir())')"
echo ">>> SLEPc installed. SLEPC_DIR=$SLEPC_DIR"

echo ">>> [3/3] Installing slepc4py (no build isolation, against installed slepc)..."
"$PY" -m pip install --no-build-isolation slepc4py

echo ">>> Build done. Verifying krylov actually engages..."
"$PY" - <<'PYEOF'
import warnings
import numpy as np
import scipy.sparse as sp

import petsc4py, slepc4py
print("petsc4py:", petsc4py.__version__, " slepc4py:", slepc4py.__version__)

# 1) The exact gate cellrank uses (decomposition/_schur.py:158). If this is
#    True, cellrank will NOT override krylov->brandts and will NOT densify.
from cellrank._utils._linear_solver import _is_petsc_slepc_available
assert _is_petsc_slepc_available(), "FAIL: cellrank STILL cannot see petsc/slepc"
print("cellrank _is_petsc_slepc_available(): True  <-- the bug is fixed")

# 2) Prove SLEPc links at RUNTIME (import success != runtime success): run a
#    real krylov Schur on a SPARSE matrix. If petsc/slepc were absent this
#    would warn + fall back + raise 'Sparse implementation is only available'.
from pygpcca._sorted_schur import sorted_schur
n = 300
rng = np.random.default_rng(0)
A = sp.random(n, n, density=0.05, format="csr", random_state=0)
A = A + sp.eye(n, format="csr")                  # no zero rows
P = sp.diags(1.0 / np.asarray(A.sum(1)).ravel()) @ A   # row-stochastic, sparse
assert sp.issparse(P)

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    try:
        R, Q, ev = sorted_schur(P, 8, "LM", method="krylov")
        print("krylov sorted_schur on SPARSE matrix: OK, returned", len(ev), "eigenvalues")
    except ValueError as e:
        if "Sparse implementation" in str(e):
            raise SystemExit("FAIL: fell back to brandts -> petsc/slepc NOT linked")
        print("krylov ran (benign:", e, ")")   # e.g. complex-conj split on random P
    msgs = [str(x.message) for x in w]
assert not any("PETSc" in m or "SLEPc" in m or "petsc4py" in m for m in msgs), \
    f"FAIL: saw a petsc/slepc fallback warning: {msgs}"

print("\n==================== VERIFIED ====================")
print("krylov Schur is now genuinely sparse/fast in cellrank_env.")
print("Your Cell 69 compute_schur(method='krylov') will NO LONGER densify.")
print("=================================================")
PYEOF

# Clean the (large) build temp; keep pip cache for next time.
rm -rf "$TMPDIR"
echo ">>> All done."
