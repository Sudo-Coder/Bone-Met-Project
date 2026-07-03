#!/usr/bin/env bash
# =============================================================================
# verify_on_hpc.sh  —  DESTINATION CLAUDE side (runs ON the HPC, in ~/prostate-cancer)
#
# Confirms the pushed data arrived intact by comparing what's on disk here
# against MANIFEST.tsv (produced by push_to_hpc.sh on the local machine).
#
#   ./hpc_transfer/verify_on_hpc.sh              # size + count check (fast)
#   ./hpc_transfer/verify_on_hpc.sh --checksum   # + sha256 of *.h5ad / *.h5 / *.pt
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # ~/prostate-cancer
MANIFEST="${ROOT}/hpc_transfer/MANIFEST.tsv"
cd "${ROOT}"

[[ -f "${MANIFEST}" ]] || { echo "ERROR: ${MANIFEST} not found — was the push completed?" >&2; exit 1; }

echo "=============================================================="
echo " Verifying ~/prostate-cancer against MANIFEST.tsv"
echo "=============================================================="

miss=0; sizebad=0; ok=0
while IFS=$'\t' read -r relpath bytes mtime; do
  [[ -z "${relpath}" ]] && continue
  if [[ ! -f "${relpath}" ]]; then
    echo "  MISSING : ${relpath}"; miss=$((miss+1)); continue
  fi
  actual=$(stat -c '%s' "${relpath}")
  if [[ "${actual}" != "${bytes}" ]]; then
    echo "  SIZE≠   : ${relpath}  (manifest ${bytes} / here ${actual})"
    sizebad=$((sizebad+1)); continue
  fi
  ok=$((ok+1))
done < "${MANIFEST}"

exp=$(wc -l < "${MANIFEST}")
echo "--------------------------------------------------------------"
printf " expected %s files | ok %s | missing %s | size-mismatch %s\n" \
       "${exp}" "${ok}" "${miss}" "${sizebad}"

# Flag anything present here but NOT in the manifest (excluding our own files)
extra=$(comm -13 \
        <(cut -f1 "${MANIFEST}" | sort) \
        <(find . -type f ! -path './hpc_transfer/*' -printf '%P\n' | sort) | wc -l)
[[ "${extra}" -gt 0 ]] && echo " note: ${extra} extra file(s) on HPC not in manifest (ok if pre-existing)"

# ----- optional deep checksum of the heavy analysis inputs ------------------
if [[ "${1:-}" == "--checksum" ]]; then
  echo "--------------------------------------------------------------"
  echo " sha256 of large inputs (*.h5ad *.h5 *.pt):"
  find . -type f \( -name '*.h5ad' -o -name '*.h5' -o -name '*.pt' \) \
       ! -path './hpc_transfer/*' -print0 |
    while IFS= read -r -d '' f; do
      printf '  %s  %s\n' "$(sha256sum "$f" | cut -c1-16)…" "${f#./}"
    done
  echo " (re-run the same command locally and compare the digests)"
fi

echo "=============================================================="
if [[ "${miss}" -eq 0 && "${sizebad}" -eq 0 ]]; then
  echo " ✅ VERIFIED — all ${exp} files present at expected sizes."
  echo "    Data is ready in ~/prostate-cancer. Safe to run analyses."
  exit 0
else
  echo " ❌ INCOMPLETE — re-run push_to_hpc.sh on the local machine;"
  echo "    rsync will resend only the missing/short files."
  exit 1
fi
