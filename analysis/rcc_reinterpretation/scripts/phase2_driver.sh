#!/usr/bin/env bash
set -u
cd /autofs/projects-t3/hussain/scProj
export MAMBA_ROOT_PREFIX=/autofs/projects-t3/hussain/scProj/.mamba TMPDIR=/autofs/scratch/arifai/tmp
LOG=analysis/rcc_reinterpretation/outputs/logs
TAB=analysis/rcc_reinterpretation/outputs/tables
RUN="bin/micromamba run -p envs/r_cellchat"
VENV=envs/rcc_reinterp_venv/bin/python
mark(){ echo "$(date '+%F %T') $*"; }

mark "waiting for CellChat/nichenetr/copykat to install..." > $LOG/phase2_driver.log
for i in $(seq 1 240); do   # up to ~2h
  if $RUN Rscript -e 'q(status=if(all(sapply(c("CellChat","nichenetr","copykat"),requireNamespace,quietly=TRUE)))0 else 1)' >/dev/null 2>&1; then
    mark "R packages ready" >> $LOG/phase2_driver.log; break
  fi
  sleep 30
done
if ! $RUN Rscript -e 'q(status=if(all(sapply(c("CellChat","nichenetr","copykat"),requireNamespace,quietly=TRUE)))0 else 1)' >/dev/null 2>&1; then
  mark "FAILED: R packages never installed" >> $LOG/phase2_driver.log; echo FAILED_PACKAGES > $TAB/PHASE2_STATUS; exit 1
fi

mark "== 20b CopyKAT ==" >> $LOG/phase2_driver.log
[ -f $TAB/copykat_predictions.csv ] || $RUN Rscript analysis/rcc_reinterpretation/scripts/20b_copykat.R >> $LOG/phase2_driver.log 2>&1
mark "== 13 CellChat ==" >> $LOG/phase2_driver.log
[ -f $TAB/cellchat_tam_LR_tumor.csv ] || $RUN Rscript analysis/rcc_reinterpretation/scripts/13_cellchat.R >> $LOG/phase2_driver.log 2>&1
mark "== 14 NicheNet ==" >> $LOG/phase2_driver.log
[ -f $TAB/nichenet_prespecified_ranks.csv ] || $RUN Rscript analysis/rcc_reinterpretation/scripts/14_nichenet.R >> $LOG/phase2_driver.log 2>&1
mark "== 15 axis support ==" >> $LOG/phase2_driver.log
$VENV analysis/rcc_reinterpretation/scripts/15_axis_support.py >> $LOG/phase2_driver.log 2>&1

S=OK
for f in cellchat_tam_LR_tumor.csv nichenet_prespecified_ranks.csv axis_support_table.csv; do
  [ -f $TAB/$f ] || S="PARTIAL(missing $f)"
done
[ -f $TAB/copykat_predictions.csv ] || S="$S;copykat_missing"
mark "DONE status=$S" >> $LOG/phase2_driver.log
echo "$S" > $TAB/PHASE2_STATUS
