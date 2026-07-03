#!/usr/bin/env bash
# =============================================================================
# push_to_hpc.sh  —  LOCAL CLAUDE side of the Bone-Met-Project -> HPC transfer
#
# Pushes the contents of Bone-Met-Project into ~/prostate-cancer on the HPC
# (mileena) using rsync over SSH. Resumable, restartable, and safe to re-run:
# rsync only sends what's missing or changed. A manifest is generated and
# shipped alongside the data so the DESTINATION CLAUDE can verify receipt.
#
# Usage:
#   ./push_to_hpc.sh <REMOTE_USER>                 # e.g. ./push_to_hpc.sh jdoe
#   REMOTE_USER=jdoe ./push_to_hpc.sh
#   ./push_to_hpc.sh jdoe --dry-run               # preview, transfer nothing
#   ./push_to_hpc.sh jdoe --checksum             # slower, verifies by checksum
# =============================================================================
set -euo pipefail

# ----- Config (override via env) --------------------------------------------
REMOTE_USER="${REMOTE_USER:-${1:-}}"
REMOTE_HOST="${REMOTE_HOST:-mileena.igs.umaryland.edu}"
REMOTE_DIR="${REMOTE_DIR:-prostate-cancer}"     # lands in ~/prostate-cancer
SRC_DIR="${SRC_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MANIFEST="${SRC_DIR}/hpc_transfer/MANIFEST.tsv"

# Strip a recognised flag if it landed in $1 instead of the username
[[ "${REMOTE_USER}" == --* ]] && REMOTE_USER=""

if [[ -z "${REMOTE_USER}" ]]; then
  echo "ERROR: no REMOTE_USER given." >&2
  echo "  usage: $0 <remote_user>   (or: REMOTE_USER=you $0)" >&2
  exit 2
fi

# Parse optional flags (any position)
DRY_RUN=""; CHECKSUM=""
for a in "$@"; do
  case "$a" in
    --dry-run)  DRY_RUN="--dry-run" ;;
    --checksum) CHECKSUM="--checksum" ;;
  esac
done

REMOTE="${REMOTE_USER}@${REMOTE_HOST}"

# ----- Exclusions -----------------------------------------------------------
# Per decision: skip ONLY the leftover partial-write temp files.
EXCLUDES=(
  --exclude='.integrated.h5ad.*'   # 4.4G interrupted-write leftover
  --exclude='*.tmp'
  --exclude='hpc_transfer/MANIFEST.remote.tsv'   # produced on the HPC side
)

echo "=============================================================="
echo " Bone-Met-Project  ->  ${REMOTE}:~/${REMOTE_DIR}/"
echo " Source : ${SRC_DIR}/"
echo " Excl.  : .integrated.h5ad.*  *.tmp"
[[ -n "$DRY_RUN"  ]] && echo " Mode   : DRY-RUN (no data will be sent)"
[[ -n "$CHECKSUM" ]] && echo " Mode   : CHECKSUM verify (slower)"
echo "=============================================================="

# ----- 1. Pre-flight: reach the host, make the target dir -------------------
echo "[1/4] Checking SSH to ${REMOTE} ..."
if ! ssh -o ConnectTimeout=15 "${REMOTE}" "mkdir -p ~/${REMOTE_DIR}/hpc_transfer && echo connected:\$(hostname)"; then
  echo "ERROR: cannot SSH to ${REMOTE}." >&2
  echo "  - Check the username is correct." >&2
  echo "  - If it asks for a password, that's fine when you run this by hand." >&2
  echo "  - To use key auth, install your key:  ssh-copy-id ${REMOTE}" >&2
  exit 1
fi

# ----- 2. Build the local manifest (relative path, bytes, mtime) ------------
echo "[2/4] Building local manifest -> ${MANIFEST}"
( cd "${SRC_DIR}" && \
  find . -type f \
      ! -name '.integrated.h5ad.*' ! -name '*.tmp' \
      ! -path './hpc_transfer/MANIFEST*.tsv' \
      -printf '%P\t%s\t%T@\n' | sort ) > "${MANIFEST}"
LOCAL_N=$(wc -l < "${MANIFEST}")
LOCAL_BYTES=$(awk -F'\t' '{s+=$2} END{print s}' "${MANIFEST}")
printf '   %s files, %s bytes (%.1f GiB)\n' "${LOCAL_N}" "${LOCAL_BYTES}" \
       "$(awk "BEGIN{print ${LOCAL_BYTES}/1073741824}")"

# ----- 3. The transfer ------------------------------------------------------
echo "[3/4] rsync ..."
# -a archive, -h human, --partial keep partial files for resume,
# --append-verify resume big files by checksumming the existing head,
# --info=progress2 single overall progress bar, --stats summary.
rsync -a -h --partial --append-verify --info=progress2 --stats \
      ${DRY_RUN} ${CHECKSUM} \
      "${EXCLUDES[@]}" \
      -e "ssh -o ConnectTimeout=15" \
      "${SRC_DIR}/" \
      "${REMOTE}:~/${REMOTE_DIR}/"

# ----- 4. Ship the manifest so the destination Claude can verify ------------
if [[ -z "$DRY_RUN" ]]; then
  echo "[4/4] Sending manifest + verify script ..."
  rsync -a "${MANIFEST}" \
        "${SRC_DIR}/hpc_transfer/verify_on_hpc.sh" \
        "${SRC_DIR}/hpc_transfer/DESTINATION_CLAUDE_HANDOFF.md" \
        "${REMOTE}:~/${REMOTE_DIR}/hpc_transfer/"
  ssh "${REMOTE}" "chmod +x ~/${REMOTE_DIR}/hpc_transfer/verify_on_hpc.sh"
else
  echo "[4/4] (dry-run) manifest not sent."
fi

echo
echo "=============================================================="
echo " DONE (local side)."
[[ -z "$DRY_RUN" ]] && cat <<EOF
 Next: tell the DESTINATION CLAUDE on ${REMOTE_HOST} to run:
     cd ~/${REMOTE_DIR} && ./hpc_transfer/verify_on_hpc.sh
 or just paste it the file:
     ~/${REMOTE_DIR}/hpc_transfer/DESTINATION_CLAUDE_HANDOFF.md
EOF
echo "=============================================================="
