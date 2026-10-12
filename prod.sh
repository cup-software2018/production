#!/bin/bash
# Produce one subrun of one detector, here and now.
#
#   ./prod.sh <fadc|iadc> <run> <subrun> [boards]
#
#   ./prod.sh fadc 263 1            FADC1 and FADC2 subrun 1
#   ./prod.sh fadc 263 1 2          FADC2 only
#   ./prod.sh iadc 263 1            IADC subrun 1 (boards are ignored)
#
# This is where the paths are decided; produce_fadc / produce_iadc only read
# the raw file they are given and write the output file they are given.
# Environment, with defaults:
#
#   RAW  /data/amore2/data/RAW     raw data, <RAW>/<run:06d>/<stream>_<run:06d>.root.<sub:05d>
#   OUT  <here>/out                output; logs go to <OUT>/log
#   CFG  <here>/prod_<det>.yaml    settings
#   BIN  <here>/produce_<det>, or <here>/build/produce_<det> if that is absent
#
# Output <OUT>/prod_<stream>_<run:06d>.<sub:05d>.root and log
# <OUT>/log/<stream>_<sub:05d>.log, for stream FADC<b> or IADC.  A stream with
# no raw file for this subrun is skipped with a note.  Exits 1 if any stream
# fails, 3 if none has input.
#
# prod_all.sh runs this for every detector of a subrun, and submit.sh runs
# prod_all.sh on the cluster, one Slurm task per subrun.  Run by hand,
# build first (cmake --build build), and mind that out/ here links to
# /data/cupsoft, which this machine cannot write to: give OUT=... instead.

set -u

DET=${1:-}
RUN=${2:-}
SUB=${3:-}
case "$DET" in fadc|iadc) ;; *) DET= ;; esac
if [ -z "$DET" ] || ! [[ "$RUN" =~ ^[0-9]+$ && "$SUB" =~ ^[0-9]+$ ]]; then
  echo "usage: $0 <fadc|iadc> <run> <subrun> [boards]"; exit 2
fi
BOARDS=${4:-"1 2"}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RAW=${RAW:-/data/amore2/data/RAW}
OUT=${OUT:-$HERE/out}
CFG=${CFG:-$HERE/prod_$DET.yaml}
if [ -z "${BIN:-}" ]; then
  BIN=$HERE/produce_$DET
  [ -x "$BIN" ] || BIN=$HERE/build/produce_$DET
fi
RUNSTR=$(printf "%06d" "$((10#$RUN))")
SUBSTR=$(printf "%05d" "$((10#$SUB))")

[ -x "$BIN" ] || { echo "no $BIN; build it with cmake first"; exit 1; }
[ -f "$CFG" ] || { echo "no config $CFG"; exit 1; }
mkdir -p "$OUT/log" || exit 1

# thisroot.sh and setup_cupdm.sh append to these, so they must exist under `set -u`.
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH:-}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:-}

# The program's RPATH finds the libraries; the libRawObjs dictionary still
# wants its headers on ROOT_INCLUDE_PATH, or cling prints "Missing FileEntry"
# for every class it autoloads.
source /home/cupsoft/Products/ROOT/V6.32.20_3.11.14/bin/thisroot.sh
source /home/cupsoft/Products/CUPDATAMODEL/V1.1/setup_cupdm.sh >/dev/null

if [ "$DET" = iadc ]; then
  streams=IADC
else
  streams=$(for b in $BOARDS; do echo "FADC$b"; done)
fi

status=0
ran=0
for s in $streams; do
  in=$RAW/$RUNSTR/${s}_$RUNSTR.root.$SUBSTR
  out=$OUT/prod_${s}_$RUNSTR.$SUBSTR.root
  log=$OUT/log/${s}_$SUBSTR.log
  if [ ! -e "$in" ]; then
    echo "note: no $in, skipped"
    continue
  fi
  {
    echo "host $(hostname)  ${SLURM_ARRAY_JOB_ID:+job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}  }$(date '+%F %T')"
    "$BIN" -c "$CFG" -o "$out" "$in"
    rc=$?
    echo "exit $rc  $(date '+%F %T')"
  } >"$log" 2>&1
  echo "$s run $RUN subrun $SUB: exit $rc, log $log"
  [ "$rc" -eq 0 ] || status=1
  ran=$((ran + 1))
done
[ "$ran" -gt 0 ] || { echo "no raw input for $DET run $RUN subrun $SUB in $RAW/$RUNSTR"; exit 3; }
exit $status
