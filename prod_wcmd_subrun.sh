#!/bin/bash
# Produce one subrun of the WCMD (FADC), here and now.
#
#   ./prod_wcmd_subrun.sh <subrun> -r <run> [-o outdir] [-c config] [boards]
#
#   ./prod_wcmd_subrun.sh 1 -r 263              FADC1 and FADC2 subrun 1, into this directory
#   ./prod_wcmd_subrun.sh 1 -r 263 -o /tmp/t    the same, into /tmp/t
#   ./prod_wcmd_subrun.sh 1 -r 263 2            FADC2 only
#
#   -r run      run number (required)
#   -o outdir   output directory, default the directory of this script (<here>)
#   -c config   settings, default <here>/prod_fadc.yaml
#
# This is where the paths are decided; produce_fadc only reads the raw file
# it is given and writes the output file it is given.  Environment, with
# defaults:
#
#   RAW  /data/amore2/data/RAW     raw data, <RAW>/<run:06d>/FADC<b>_<run:06d>.root.<sub:05d>
#   BIN  <here>/produce_fadc, or <here>/build/produce_fadc if that is absent
#   CONF /data/amore2/data/CONFIG  run configs, <CONF>/<run:06d>.yml
#
# The windows follow the run's FADC DLY, read from the run config, and reach
# prod_fadc.yaml as ${FADC_BASE_HI}, ${FADC_INT_LO}, ${FADC_INT_HI}.
#
# Output <outdir>/prod_FADC<b>_<run:06d>.<sub:05d>.root and log
# <outdir>/log/FADC<b>_<sub:05d>.log.  A board with no raw file for this
# subrun is skipped with a note.  Exits 1 if any board fails, 3 if none has
# input.  Build first (cmake --build build).

set -u

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
usage() { echo "usage: $0 <subrun> -r <run> [-o outdir] [-c config] [boards]"; exit 2; }
RUN=
OUT=$HERE
CFG=$HERE/prod_fadc.yaml
args=()
while [ $# -gt 0 ]; do
  case "$1" in
    -r|-o|-c)
      [ $# -ge 2 ] || usage
      case "$1" in -r) RUN=$2 ;; -o) OUT=$2 ;; -c) CFG=$2 ;; esac
      shift 2 ;;
    -*) usage ;;
    *) args+=("$1"); shift ;;
  esac
done
SUB=${args[0]:-}
BOARDS=${args[1]:-"1 2"}
[[ "$RUN" =~ ^[0-9]+$ && "$SUB" =~ ^[0-9]+$ && ${#args[@]} -le 2 ]] || usage

RAW=${RAW:-/data/amore2/data/RAW}
CONF=${CONF:-/data/amore2/data/CONFIG}
if [ -z "${BIN:-}" ]; then
  BIN=$HERE/produce_fadc
  [ -x "$BIN" ] || BIN=$HERE/build/produce_fadc
fi
RUNSTR=$(printf "%06d" "$((10#$RUN))")
SUBSTR=$(printf "%05d" "$((10#$SUB))")

[ -x "$BIN" ] || { echo "no $BIN; build it with cmake first"; exit 1; }
[ -f "$CFG" ] || { echo "no config $CFG"; exit 1; }
mkdir -p "$OUT/log" || exit 1

# Windows from the DLY, the same on every board (windows.sh).
. "$HERE/windows.sh"
fadc_windows "$CONF/$RUNSTR.yml" || exit 1
dlynote=$WINDOWS_NOTE
echo "$dlynote"

# thisroot.sh and setup_cupdm.sh append to these, so they must exist under `set -u`.
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH:-}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:-}

# The program's RPATH finds the libraries; the libRawObjs dictionary still
# wants its headers on ROOT_INCLUDE_PATH, or cling prints "Missing FileEntry"
# for every class it autoloads.
source /home/cupsoft/Products/ROOT/V6.32.20_3.11.14/bin/thisroot.sh
source /home/cupsoft/Products/CUPDATAMODEL/V1.1/setup_cupdm.sh >/dev/null

status=0
ran=0
for b in $BOARDS; do
  s=FADC$b
  in=$RAW/$RUNSTR/${s}_$RUNSTR.root.$SUBSTR
  out=$OUT/prod_${s}_$RUNSTR.$SUBSTR.root
  log=$OUT/log/${s}_$SUBSTR.log
  if [ ! -e "$in" ]; then
    echo "note: no $in, skipped"
    continue
  fi
  {
    echo "host $(hostname)  $(date '+%F %T')"
    echo "$dlynote"
    "$BIN" -c "$CFG" -o "$out" "$in"
    rc=$?
    echo "exit $rc  $(date '+%F %T')"
  } >"$log" 2>&1
  echo "$s run $RUN subrun $SUB: exit $rc, log $log"
  [ "$rc" -eq 0 ] || status=1
  ran=$((ran + 1))
done
[ "$ran" -gt 0 ] || { echo "no raw FADC input for run $RUN subrun $SUB in $RAW/$RUNSTR"; exit 3; }
exit $status
