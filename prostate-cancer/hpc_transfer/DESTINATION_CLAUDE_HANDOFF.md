# Handoff: Local Claude → Destination Claude (mileena HPC)

You are the **destination Claude** running on `mileena.igs.umaryland.edu`. The
**local Claude** (on the user's WSL machine) is pushing the `Bone-Met-Project`
data to you via `rsync` over SSH. It lands in **`~/prostate-cancer/`**.

This is a **single-cell prostate-cancer bone-metastasis** analysis project
(scRNA-seq: scVI integration, CellRank, CellOracle GRNs). The heavy inputs are
several ~7 GB `.h5ad` files plus `Data/` count matrices.

## Your job

1. **Confirm the landing dir exists** (rsync/local Claude create it):
   ```bash
   ls -la ~/prostate-cancer/
   ```

2. **Wait until the push is done.** The local side reports when rsync finishes.
   Contents arrive incrementally; don't verify mid-transfer.

3. **Verify integrity** against the manifest the local Claude shipped:
   ```bash
   cd ~/prostate-cancer && ./hpc_transfer/verify_on_hpc.sh
   # deeper (checksums the big .h5ad/.h5/.pt inputs):
   ./hpc_transfer/verify_on_hpc.sh --checksum
   ```
   - ✅ "VERIFIED" → data is complete, ready to run.
   - ❌ "INCOMPLETE" → tell the user; local Claude re-runs `push_to_hpc.sh` and
     rsync resends only what's missing/short. Then verify again.

## What was intentionally NOT sent

- `.integrated.h5ad.GqjYY1` and `*.tmp` — leftover partial-write junk. Do **not**
  expect these; their absence is correct.
- Everything else in `Bone-Met-Project` (including `.git`) **is** sent.

## Notes for running analyses here

- Primary notebook per the project README: `single_cell_analysis_complete_class.ipynb`.
- Env recreation notes are in the repo: `RECREATE_ENVIRONMENTS_CLAUDE_PROMPT.md`
  and `celloracle_env_setup_notes.md`. Set up the conda/venv envs **on the HPC**
  from those before running — the Python environments were **not** transferred,
  only project files and data.
- If you need to send results back to the user's machine, mirror in reverse:
  `rsync -avhP ~/prostate-cancer/Outputs/ <localuser>@<localhost>:.../Outputs/`
  (or just leave them and the local Claude will pull).

## Coordination protocol (optional)

To signal status back, write a line to `~/prostate-cancer/hpc_transfer/STATUS`:
```bash
echo "$(hostname): verified OK, ready to run" >> ~/prostate-cancer/hpc_transfer/STATUS
```
The local Claude can `ssh` in and read that file to know you're done.
