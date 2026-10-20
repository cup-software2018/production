#!/bin/bash
# PSMD gain matching of one run, start to end.
#
#   ./run_psmd_gain.sh <run> [prod]        prod: production directory, default /data/amore2/data/PROD
#
# Runs 1_ped_hist, 2_dark_gain, 3_muon_cal and the three checks in a scratch
# directory under <prod>/psmd_gain, keeps only
#
#   <prod>/psmd_gain/psmd_gain_<run>.txt              the table analyses read: pid g c use
#   <prod>/psmd_gain/check_dark_gain_<run>.png        g of every channel
#   <prod>/psmd_gain/check_dark_gain_1_<run>.png      dark-pulse spectra of PID 1 with the fitted p.e.
#   <prod>/psmd_gain/check_gain_match_<run>.png       clean-muon response before / after g
#   <prod>/psmd_gain/check_muon_cal_<run>.png         c of every channel
#   <prod>/psmd_gain/check_odd_match_<run>.png        matching per channel on the odd subruns, c from the even
#   <prod>/psmd_gain/psmd_gain_<run>.log              the output of every step
#
# and removes the intermediate files (histograms, npz).  Exits non-zero, keeping
# the scratch directory, if a step fails.

set -u
RUN=${1:-}
PROD=${2:-/data/amore2/data/PROD}
if ! [[ "$RUN" =~ ^[0-9]+$ ]]; then
  echo "usage: $0 <run> [prod]"; exit 2
fi
RUN=$((10#$RUN))
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROD=$(cd "$PROD" && pwd) || exit 1
OUT=$PROD/psmd_gain
mkdir -p "$OUT" || exit 1
WORK=$(mktemp -d "$OUT/work_${RUN}_XXXXXX") || exit 1
LOG=$OUT/psmd_gain_${RUN}.log
PY=${PY:-python}

step() {
  echo "== $*" | tee -a "$LOG"
  if ! (cd "$WORK" && "$PY" "$HERE/$1" "${@:2}") >>"$LOG" 2>&1; then
    echo "FAILED: $*  (scratch kept: $WORK, log: $LOG)"
    exit 1
  fi
}

: >"$LOG"
echo "PSMD gain matching, run $RUN, production $PROD  $(date '+%F %T')" | tee -a "$LOG"
step 1_ped_hist.py "$RUN" --prod "$PROD" --dir "$WORK"
step 2_dark_gain.py "$RUN" --dir "$WORK"
step 3_muon_cal.py "$RUN" --prod "$PROD" --dir "$WORK"
step check_dark_gain.py "$RUN" --dir "$WORK"
step check_muon_cal.py "$RUN" --dir "$WORK"
step check_gain_match.py hist "$RUN" --prod "$PROD" --dir "$WORK"
step check_gain_match.py plot "$RUN" --dir "$WORK"

mv "$WORK/psmd_gain_${RUN}.txt" "$WORK"/check_*_"${RUN}".png "$OUT"/ || exit 1
rm -rf "$WORK"
echo "done $(date '+%F %T'): $OUT/psmd_gain_${RUN}.txt and check_*_${RUN}.png" | tee -a "$LOG"
