#!/bin/bash
# WCMD gain matching of one run, start to end.
#
#   ./run_wcmd_gain.sh <run> [prod] [FIRST:LAST]   prod: production directory, default /data/amore2/data/PROD
#                                                  FIRST:LAST: only these subruns (for tests)
#
# Runs 1_dark_hist, 2_spe_fit and the check in a scratch directory under <prod>/wcmd_gain, keeps only
#
#   <prod>/wcmd_gain/wcmd_gain_<run>.txt          the table analyses read: pid g ... use
#   <prod>/wcmd_gain/check_spe_fit_<run>.png      the 47 SPE fits
#   <prod>/wcmd_gain/check_gain_<run>.png         g of every tank PMT, fit quality
#   <prod>/wcmd_gain/wcmd_gain_<run>.log          the output of every step
#
# and removes the intermediate files (histograms, npz).  Exits non-zero, keeping the scratch
# directory, if a step fails.

set -u
RUN=${1:-}
PROD=${2:-/data/amore2/data/PROD}
SUBS=${3:-}
if ! [[ "$RUN" =~ ^[0-9]+$ ]]; then
  echo "usage: $0 <run> [prod] [FIRST:LAST]"; exit 2
fi
RUN=$((10#$RUN))
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROD=$(cd "$PROD" && pwd) || exit 1
OUT=$PROD/wcmd_gain
mkdir -p "$OUT" || exit 1
WORK=$(mktemp -d "$OUT/work_${RUN}_XXXXXX") || exit 1
LOG=$OUT/wcmd_gain_${RUN}.log
PY=${PY:-python}

step() {
  echo "== $*" | tee -a "$LOG"
  if ! (cd "$WORK" && "$PY" "$HERE/$1" "${@:2}") >>"$LOG" 2>&1; then
    echo "FAILED: $*  (scratch kept: $WORK, log: $LOG)"
    exit 1
  fi
}

: >"$LOG"
echo "WCMD gain matching, run $RUN, production $PROD ${SUBS:+subruns $SUBS} $(date '+%F %T')" | tee -a "$LOG"
step 1_dark_hist.py "$RUN" --prod "$PROD" --dir "$WORK" ${SUBS:+--subruns "$SUBS"}
step 2_spe_fit.py "$RUN" --dir "$WORK"
step check_spe_fit.py "$RUN" --dir "$WORK"

mv "$WORK/wcmd_gain_${RUN}.txt" "$WORK"/check_*_"${RUN}".png "$OUT"/ || exit 1
rm -rf "$WORK"
echo "done $(date '+%F %T'): $OUT/wcmd_gain_${RUN}.txt and check_*_${RUN}.png" | tee -a "$LOG"
