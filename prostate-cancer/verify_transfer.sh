#!/usr/bin/env bash
# verify_transfer.sh — completeness check for the Bone-Met-Project -> prostate-cancer transfer.
#
# Run the SAME script on BOTH ends (identical output = complete, byte-for-byte transfer):
#   LAPTOP (source):   ./verify_transfer.sh /path/to/Bone-Met-Project
#   HPC (destination): ./verify_transfer.sh /autofs/projects-t3/hussain/prostate-cancer
#
# It prints a per-file count, a total byte count, and a single deterministic
# sha256 "manifest digest" over (relative-path + size + content-hash) of every
# regular file. If the digest matches on both ends, the transfer is verified.
set -euo pipefail

ROOT="${1:?usage: verify_transfer.sh <dir>}"
ROOT="${ROOT%/}"
[ -d "$ROOT" ] || { echo "ERROR: not a directory: $ROOT" >&2; exit 1; }

cd "$ROOT"

# Deterministic file list: all regular files, path-sorted, excluding the
# script itself and common noise so both ends agree.
mapfile -d '' FILES < <(
  find . -type f \
    ! -name 'verify_transfer.sh' \
    ! -name 'PROMPT_FOR_LOCAL_CLAUDE.txt' \
    ! -path './.git/*' \
    ! -name '.DS_Store' \
    -print0 | sort -z
)

COUNT=${#FILES[@]}
if [ "$COUNT" -eq 0 ]; then
  echo "root:            $ROOT"
  echo "files:           0   (empty — nothing transferred yet)"
  exit 0
fi

BYTES=$(printf '%s\0' "${FILES[@]}" | du -cb --files0-from=- 2>/dev/null | tail -1 | cut -f1)

# Manifest digest: for each file emit "sha256  size  ./relative/path", then
# hash that whole listing. Independent of timestamps, owners, or transfer order.
DIGEST=$(
  for f in "${FILES[@]}"; do
    h=$(sha256sum "$f" | cut -d' ' -f1)
    s=$(stat -c%s "$f")
    printf '%s  %s  %s\n' "$h" "$s" "$f"
  done | sha256sum | cut -d' ' -f1
)

echo "root:            $ROOT"
echo "files:           $COUNT"
echo "total bytes:     $BYTES"
echo "manifest digest: $DIGEST"
