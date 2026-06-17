# CellOracle Environment — Reproducibility Recipe

**Created:** 2026-06-03
**System:** WSL2 Ubuntu 22.04.5 LTS, x86-64
**Conda:** Miniconda (conda 26.3.2), channel = conda-forge only, strict priority
**Env:** `celloracle_env`, Python 3.10

This is the from-scratch recipe for rebuilding the CellOracle environment, including
every fix needed beyond a plain `pip install celloracle`. CellOracle is fragile because
several of its dependencies (velocyto, gimmemotifs, pybedtools) compile legacy C/C++ from
source and were written against older toolchains.

---

## Verified working versions (end state)

| Package      | Version  |
|--------------|----------|
| python       | 3.10.20  |
| celloracle   | 0.20.0   |
| numpy        | 1.26.4   |
| pandas       | 1.5.3    |  ← CellOracle pins `pandas<=1.5.3`
| scanpy       | 1.11.0   |
| anndata      | 0.10.8   |
| numba        | 0.65.1   |
| scikit-learn | 1.5.2    |
| ipykernel    | 7.2.0    |

---

## Full rebuild — copy/paste

```bash
# 0. Channel setup (one-time, avoids Anaconda ToS; conda-forge is free)
conda config --add channels conda-forge
conda config --set channel_priority strict

# 1. Create env + pip (conda-forge python does NOT bundle pip)
conda create -n celloracle_env python=3.10 -y -c conda-forge --override-channels
conda activate celloracle_env
conda install -n celloracle_env pip -y -c conda-forge --override-channels

# 2. Build toolchain + zlib headers (no system gcc on a fresh WSL box; sudo-free)
conda install -n celloracle_env c-compiler cxx-compiler zlib -y -c conda-forge --override-channels
conda activate celloracle_env   # re-activate so CC/CXX point at the conda toolchain

# 3. Pre-seed numpy<2 + Cython, then build velocyto WITHOUT build isolation
#    (velocyto's setup.py does `import numpy` at build time; PEP517 isolation hides it)
pip install "numpy<2" cython
pip install velocyto --no-build-isolation

# 4. Install CellOracle. Relax GCC-14 hard-errors for gimmemotifs' legacy C.
#    (GCC 14 promoted implicit-int / int-conversion from warnings to errors)
export CFLAGS="$CFLAGS -Wno-error=implicit-int -Wno-error=int-conversion -Wno-error=implicit-function-declaration"
pip install celloracle ipykernel

# 5. Provide pkg_resources for gimmemotifs (conda-forge python omits setuptools;
#    pin <81 because setuptools 81+ removes pkg_resources)
pip install "setuptools<81"

# 6. Verify
python -c "import celloracle; print('celloracle', celloracle.__version__)"
```

---

## The five fixes (why each was needed)

A plain `pip install celloracle` failed repeatedly on a clean conda-forge env. Each fix
addressed a distinct, separate failure — none were CellOracle version conflicts; the pip
resolver itself succeeded. They were all **build-tooling / packaging** problems:

1. **pip missing.** conda-forge's Python build does not include pip by default.
   → `conda install pip`.

2. **velocyto build failed: `ModuleNotFoundError: No module named 'numpy'`.**
   velocyto's old-style `setup.py` imports numpy at build time, but pip's isolated
   PEP-517 build environment doesn't have numpy yet.
   → Pre-install `numpy<2` + `cython`, then `pip install velocyto --no-build-isolation`
   so the build sees the env's numpy. (`numpy<2` because numba / pandas 1.5.3 require
   the numpy 1.x ABI — numpy 2.0 breaks the CellOracle stack.)

3. **No C compiler: `command 'gcc' failed: No such file or directory`.**
   Fresh WSL Ubuntu has no compiler, and sudo needed a password (couldn't apt-install
   build-essential). Solved sudo-free via conda.
   → `conda install c-compiler` (conda-forge GCC toolchain, sets `$CC` on activation).

4. **No C++ compiler + missing zlib headers.**
   pybedtools needs g++ (`command 'g++' failed`) and zlib (`fatal error: zlib.h: No such
   file or directory`).
   → `conda install cxx-compiler zlib`. The conda compiler auto-adds
   `-I$CONDA_PREFIX/include`, so it finds `zlib.h` once zlib is installed.

5. **gimmemotifs legacy C rejected by GCC 14, then `pkg_resources` missing at import.**
   - Build: `error: type of 'len1' defaults to 'int' [-Wimplicit-int]` and
     `[-Wint-conversion]`. GCC 14 made these hard errors by default.
     → `export CFLAGS="... -Wno-error=implicit-int -Wno-error=int-conversion
       -Wno-error=implicit-function-declaration"` (downgrades back to warnings; behavior
       unchanged — it's pre-C99 source).
   - Import: `ModuleNotFoundError: No module named 'pkg_resources'`. gimmemotifs imports
     pkg_resources, which ships with setuptools — not bundled by conda-forge python.
     → `pip install "setuptools<81"` (81+ removes pkg_resources).

---

## Notes

- On first import, gimmemotifs prints `INFO - Using included version of ...` and
  `WARNING - <tool> not found` for MEME, Homer, DREME, etc. These are **optional**
  external de-novo motif-discovery tools. CellOracle's default motif scanning uses the
  bundled scanner and does not need them. The config is written to
  `~/.config/gimmemotifs/gimmemotifs.cfg`.
- A `pkg_resources is deprecated` UserWarning is harmless.
- Kept install logs: `~/celloracle_install*.log`, `~/velocyto_install*.log`.
- `cellrank` lives in a SEPARATE env (`cellrank_env`, Python 3.11) — do not mix; their
  dependencies conflict.
