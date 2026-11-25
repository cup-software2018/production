#!/bin/bash
# Produce one subrun of the PSMD (IADC), here and now.
#
#   ./prod_psmd_subrun.sh <subrun> -r <run> [-o outdir] [-c config]
#
#   ./prod_psmd_subrun.sh 1 -r 343              subrun 1, into this directory
#   ./prod_psmd_subrun.sh 1 -r 343 -o /tmp/t    the same, into /tmp/t
#
#   -r run      run number (required)
#   -o outdir   output directory, default the directory of this script (<here>)
#   -c config   settings, default <here>/prod_iadc.yaml
#
# This is where the paths are decided; produce_iadc only reads the raw file
# it is given and writes the output file it is given.  Environment, with
# defaults:
#
#   RAW  /data/amore2/data/RAW     raw data, <RAW>/<run:06d>/IADC_<run:06d>.root.<sub:05d>
#   BIN  <here>/produce_iadc, or <here>/build/produce_iadc if that is absent
#   CONF /data/amore2/data/CONFIG  run configs, <CONF>/<run:06d>.yml
#   IADC_THREADS 1                 threads of produce_iadc (-j)
#
# The windows follow the run's IADC DLY, read from the run config, and reach
# prod_iadc.yaml as ${IADC_BASE_HI}, ${IADC_INT_LO}, ${IADC_INT_HI}.
#
# Output <outdir>/prod_IADC_<run:06d>.<sub:05d>.root and log
# <outdir>/log/IADC_<sub:05d>.log.  Exits 1 if it fails, 3 if there is no raw
# file for this subrun.  Build first (cmake --build build).

set -u

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
usage() { echo "usage: $0 <subrun> -r <run> [-o outdir] [-c config]"; exit 2; }
RUN=
OUT=$HERE
CFG=$HERE/prod_iadc.yaml
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
[[ "$RUN" =~ ^[0-9]+$ && "$SUB" =~ ^[0-9]+$ && ${#args[@]} -le 1 ]] || usage

RAW=${RAW:-/data/amore2/data/RAW}
CONF=${CONF:-/data/amore2/data/CONFIG}
if [ -z "${BIN:-}" ]; then
  BIN=$HERE/produce_iadc
  [ -x "$BIN" ] || BIN=$HERE/build/produce_iadc
fi
RUNSTR=$(printf "%06d" "$((10#$RUN))")
SUBSTR=$(printf "%05d" "$((10#$SUB))")

in=$RAW/$RUNSTR/IADC_$RUNSTR.root.$SUBSTR
out=$OUT/prod_IADC_$RUNSTR.$SUBSTR.root
log=$OUT/log/IADC_$SUBSTR.log
[ -e "$in" ] || { echo "no raw IADC input for run $RUN subrun $SUB: $in"; exit 3; }
[ -x "$BIN" ] || { echo "no $BIN; build it with cmake first"; exit 1; }
[ -f "$CFG" ] || { echo "no config $CFG"; exit 1; }
mkdir -p "$OUT/log" || exit 1

# Windows from the DLY (windows.sh).
. "$HERE/windows.sh"
iadc_windows "$CONF/$RUNSTR.yml" || exit 1
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

{
  echo "host $(hostname)  $(date '+%F %T')"
  echo "$dlynote"
  "$BIN" -c "$CFG" -j "${IADC_THREADS:-1}" -o "$out" "$in"
  rc=$?
  echo "exit $rc  $(date '+%F %T')"
} >"$log" 2>&1
echo "IADC run $RUN subrun $SUB: exit $rc, log $log"
[ "$rc" -eq 0 ] || exit 1
exit 0
