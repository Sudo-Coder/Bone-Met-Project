# Task for Claude Code: Recreate 3 conda environments + Jupyter kernels

You are setting up a computational-biology Python stack on a freshly SSH'd Linux machine.
Recreate **three** conda environments exactly as specified, register a Jupyter kernel for
each, and write reproducibility lockfiles. This recipe already encodes every fix discovered
when these envs were first built — **follow the step order exactly** so each env builds in a
single clean pass. Do not improvise a different order; the ordering exists to avoid known
build failures.

## Operating rules (read first)
- Work **one logical step at a time**; after each major step, show the command output and
  confirm success before continuing. Do not chain the whole thing into one blind script.
- **Before any large download or long step, say so.** Large steps are flagged with 📦 below.
- **Do NOT install anything into `base`** — keep base for conda itself. Every package goes
  into a named env.
- **Do NOT delete or overwrite any existing user files.** Only create new conda envs,
  kernels, and the lockfiles named at the end. If a target env already exists, STOP and ask.
- Prefer **conda-forge**; use pip only where noted.
- If a step fails, **diagnose and explain before retrying** — do not blindly force-install,
  and do not use `--force-reinstall` / `--ignore-installed` without explaining why.
- If you hit a dependency *version* conflict not described here, STOP and report it with the
  exact error before proceeding.

## Pre-flight checks
Run these and confirm before starting:
```bash
uname -m                 # expect x86_64. If aarch64/arm64, change the Miniconda URL below
                         #   to the matching -aarch64 installer.
cat /etc/os-release | grep -E '^(NAME|VERSION)='   # expect a Linux distro
which conda || echo "no conda yet"                 # if conda already exists, do NOT reinstall;
                                                   #   skip Phase 0 and reuse it.
df -h ~ | tail -1        # ensure >5 GB free (scvi-tools/PyTorch alone is ~1-2 GB)
which gcc g++ || echo "no system compiler (fine — we install one via conda)"
```
This recipe assumes **Linux x86-64**. It does **not** require sudo or a system compiler —
the C/C++ toolchain is installed via conda-forge into the env that needs it.

---

## Phase 0 — Install Miniconda (skip if conda already present)

📦 Large download (~150 MB).
```bash
cd ~
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda_installer.sh
bash ~/miniconda_installer.sh -b -p ~/miniconda3        # -b = batch, only creates ~/miniconda3
~/miniconda3/bin/conda init bash                        # only this edits ~/.bashrc (clean block)
rm ~/miniconda_installer.sh
```
For all subsequent steps, activate conda non-interactively at the top of each shell:
```bash
source ~/miniconda3/etc/profile.d/conda.sh
```

### Channel config — avoids the Anaconda Terms-of-Service prompt
conda now blocks `repo.anaconda.com` default channels behind a ToS acceptance, and those
channels may require a paid license for some institutions. We sidestep it entirely by using
the free community **conda-forge** channel (standard for bioinformatics):
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```
From here on, every `conda create`/`conda install` uses `-c conda-forge --override-channels`.

Verify base works (expect conda 26.x, base python 3.13.x — base is only for conda):
```bash
conda --version
conda run -n base python --version
```

---

## Phase 1 — `celloracle_env` (Python 3.10)  ← the fragile one, do steps IN ORDER

**Why this is delicate:** CellOracle pulls dependencies (`velocyto`, `gimmemotifs`,
`pybedtools`) that compile legacy C/C++ from source and were written for older toolchains.
A naive `pip install celloracle` fails ~5 different ways. The steps below pre-empt all of
them. Target end-state versions to verify against are in the table at the bottom.

```bash
source ~/miniconda3/etc/profile.d/conda.sh

# 1.1 Create env + pip (conda-forge python does NOT bundle pip)
conda create -n celloracle_env python=3.10 -y -c conda-forge --override-channels
conda install -n celloracle_env pip -y -c conda-forge --override-channels

# 1.2 Install the FULL build toolchain + zlib headers UP FRONT (sudo-free).
#     Reason: velocyto needs gcc, pybedtools needs g++ AND zlib.h. Installing all now
#     avoids three separate failure-and-retry rounds.
conda install -n celloracle_env c-compiler cxx-compiler zlib -y -c conda-forge --override-channels

# 1.3 Activate so $CC/$CXX point at the conda toolchain (the conda compilers auto-add
#     -I$CONDA_PREFIX/include, which is how zlib.h gets found).
conda activate celloracle_env

# 1.4 Pre-seed numpy<2 + Cython, then build velocyto WITHOUT build isolation.
#     Reason: velocyto's setup.py does `import numpy` at build time, but pip's isolated
#     PEP-517 build env hides it -> ModuleNotFoundError: numpy. --no-build-isolation lets
#     the build see the env's numpy. numpy MUST be <2 (numba + pandas 1.5.3 need the 1.x ABI).
pip install "numpy<2" cython
pip install velocyto --no-build-isolation

# 1.5 Relax GCC-14 hard-errors for gimmemotifs' pre-C99 source, then install celloracle.
#     Reason: GCC 14 promoted implicit-int / int-conversion from warnings to ERRORS;
#     gimmemotifs' old C won't compile without downgrading them back to warnings.
export CFLAGS="$CFLAGS -Wno-error=implicit-int -Wno-error=int-conversion -Wno-error=implicit-function-declaration"
pip install celloracle ipykernel

# 1.6 Provide pkg_resources for gimmemotifs at import time.
#     Reason: conda-forge python omits setuptools; gimmemotifs imports pkg_resources.
#     Pin <81 because setuptools 81+ removes pkg_resources.
pip install "setuptools<81"

# 1.7 Verify (first import triggers a one-time gimmemotifs config write + INFO/WARNING
#     lines about optional motif tools like MEME/Homer — those are EXPECTED and harmless).
python -c "import celloracle; print('celloracle', celloracle.__version__)"
python -c "import numpy,pandas,scanpy,anndata,numba; print('numpy',numpy.__version__,'| pandas',pandas.__version__,'| scanpy',scanpy.__version__,'| anndata',anndata.__version__,'| numba',numba.__version__)"
```
**STOP and report** the printed versions. Expected: celloracle 0.20.0, numpy 1.26.x,
pandas 1.5.3, scanpy 1.11.x, anndata 0.10.x, python 3.10.x. If celloracle imports and the
versions match the table, Phase 1 succeeded.

---

## Phase 2 — `cellrank_env` (Python 3.11)  ← easy, modern wheels

CellRank is well-behaved (proper wheels, no exotic C builds). It uses numpy 2.x / pandas 2.x
and therefore **must stay in its own env, separate from celloracle_env** (which is pinned to
numpy 1.x / pandas 1.5.3 — they cannot coexist).

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda create -n cellrank_env python=3.11 -y -c conda-forge --override-channels
conda install -n cellrank_env pip -y -c conda-forge --override-channels
conda activate cellrank_env
pip install cellrank ipykernel          # 📦 sizable, but clean
python -c "import cellrank; print('cellrank', cellrank.__version__)"
python -c "import numpy,pandas,scanpy,anndata; print('numpy',numpy.__version__,'| pandas',pandas.__version__,'| scanpy',scanpy.__version__,'| anndata',anndata.__version__)"
```
**STOP and report.** Expected: cellrank 2.0.7, numpy 2.x, pandas 2.3.x, python 3.11.x.

---

## Phase 3 — `sc-general-analysis` (Python 3.11)  ← everyday single-cell env

General-purpose env derived from scanning two analysis notebooks. It contains the scientific
stack + scanpy/anndata/scvi-tools/gseapy. **`diffxpy` is deliberately EXCLUDED** — it is
unmaintained and depends on an old TensorFlow chain that conflicts with scvi-tools' PyTorch
stack. Do not add diffxpy here.

📦 This pulls **PyTorch (CPU, ~1-2 GB)** via scvi-tools. Flag it before running.
```bash
source ~/miniconda3/etc/profile.d/conda.sh

# 3.1 Create env with the full conda-forge stack in one solve.
#     leidenalg + python-igraph are scanpy's leiden-clustering runtime deps (not imported
#     directly but needed by sc.tl.leiden). scvi-tools provides `import scvi` and pulls CPU PyTorch.
conda create -n sc-general-analysis -y -c conda-forge --override-channels \
  python=3.11 \
  numpy pandas scipy matplotlib seaborn scikit-learn \
  jupyter ipykernel ipython \
  anndata scanpy scvi-tools leidenalg python-igraph \
  pip

# 3.2 gseapy is not reliably on conda-forge (it's bioconda) -> pip. Good PyPI wheels exist.
conda activate sc-general-analysis
pip install gseapy

# 3.3 Verify the exact imports the notebooks use (note: import name is `scvi`, not `scvi-tools`).
python - <<'PYEOF'
import importlib, sys
print("python", sys.version.split()[0])
for imp in ["numpy","pandas","scipy","matplotlib","seaborn","sklearn",
            "anndata","scanpy","scvi","gseapy","IPython"]:
    m = importlib.import_module(imp)
    print(f"{imp:12}", getattr(m,"__version__","(ok)"))
PYEOF
```
**STOP and report.** Expected: scvi 1.4.x, scanpy 1.11.x, anndata 0.12.x, numpy 2.x,
gseapy 1.2.x, python 3.11.x. (anndata/scanpy may print a harmless `__version__ deprecated`
FutureWarning.)

---

## Phase 4 — Register kernels + write lockfiles

```bash
source ~/miniconda3/etc/profile.d/conda.sh

# 4.1 One Jupyter kernel per env (names chosen to match the original machine).
conda run -n celloracle_env       python -m ipykernel install --user --name celloracle_env       --display-name "Python (celloracle_env)"
conda run -n cellrank_env         python -m ipykernel install --user --name cellrank_env         --display-name "Python (cellrank_env)"
conda run -n sc-general-analysis  python -m ipykernel install --user --name sc-general-analysis  --display-name "Python (sc-general-analysis)"
jupyter --version >/dev/null 2>&1 && jupyter kernelspec list   # or: conda run -n sc-general-analysis jupyter kernelspec list

# 4.2 Reproducibility lockfiles in ~ (one conda manifest + one pip-freeze per env).
for env in celloracle_env cellrank_env sc-general-analysis; do
  conda env export -n "$env" > ~/"${env}_env.yml"
  conda run -n "$env" python -m pip freeze > ~/"${env}_pip_freeze.txt"
done
ls -lh ~/*_env.yml ~/*_pip_freeze.txt
```

In a VS Code notebook on this machine, the kernels appear in the picker as
**Python (celloracle_env)**, **Python (cellrank_env)**, **Python (sc-general-analysis)**.
Run **Developer: Reload Window** if they don't show up immediately.

---

## Expected end-state version tables (use to confirm "exact" recreation)

**celloracle_env** — Python 3.10.20
| celloracle 0.20.0 | numpy 1.26.4 | pandas 1.5.3 | scanpy 1.11.0 | anndata 0.10.8 | numba 0.65.1 | scikit-learn 1.5.2 | ipykernel 7.2.0 |

**cellrank_env** — Python 3.11.15
| cellrank 2.0.7 | numpy 2.4.6 | pandas 2.3.3 | scanpy 1.11.5 | anndata 0.12.16 |

**sc-general-analysis** — Python 3.11.15
| scvi-tools 1.4.2 | scanpy 1.11.5 | anndata 0.12.16 | numpy 2.4.6 | pandas 2.3.3 | scipy 1.17.1 | matplotlib 3.10.9 | seaborn 0.13.2 | scikit-learn 1.9.0 | gseapy 1.2.1 |

Minor point-release drift is normal if rebuilt later. For **byte-exact** reproduction
instead of recipe-based, see below.

---

## Optional: byte-exact reproduction via lockfiles (highest fidelity)
If you copy the original machine's lockfiles to this box (`scp` the `*_env.yml` files into
`~`), you can try exact restore instead of the recipe:
```bash
conda env create -n celloracle_env -f ~/celloracle_env_env.yml
# ...repeat per env...
```
Caveats: full `conda env export` pins build strings that must exist on this platform/arch;
if `conda env create` fails to solve, fall back to the recipe phases above (which are the
robust path). The celloracle pip chain in particular is most reliably reproduced via the
Phase 1 recipe, not the yml. The `*_pip_freeze.txt` files give the exact pip versions if you
want to pin any pip install precisely.

---

## Troubleshooting quick-reference (celloracle, if a step still fails)
| Symptom | Cause | Fix |
|---|---|---|
| `No module named 'numpy'` building velocyto | numpy hidden by build isolation | ensure 1.4 ran: `pip install "numpy<2" cython` then `velocyto --no-build-isolation` |
| `command 'gcc'/'g++' failed: No such file` | no compiler | `conda install c-compiler cxx-compiler` then re-activate env |
| `fatal error: zlib.h` | missing zlib headers | `conda install zlib` (conda compiler finds it via `$CONDA_PREFIX/include`) |
| gimmemotifs `error: ... [-Wimplicit-int]` / `[-Wint-conversion]` | GCC 14 strictness | export the `CFLAGS="... -Wno-error=..."` from step 1.5 before `pip install celloracle` |
| `No module named 'pkg_resources'` on `import celloracle` | setuptools absent | `pip install "setuptools<81"` |
```
