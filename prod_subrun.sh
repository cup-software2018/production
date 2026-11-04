#!/bin/bash
# Produce one subrun of one detector: prod_wcmd_subrun.sh (fadc) or
# prod_psmd_subrun.sh (iadc), with the output in OUT.
#
#   ./prod_subrun.sh <fadc|iadc> <run> <subrun> [boards]
#
#   ./prod_subrun.sh fadc <run> 1          FADC1 and FADC2 subrun 1
#   ./prod_subrun.sh fadc <run> 1 2        FADC2 only
#   ./prod_subrun.sh iadc <run> 1          IADC subrun 1 (boards are ignored)
#
#   OUT  <here>/out                 output (-o); logs go to <OUT>/log
#   CFG  the script's default       settings (-c)
#
# RAW, BIN and CONF go through to those scripts, as do their exit codes.
# prod_all_subrun.sh runs this for every detector of a subrun, and submit.sh runs
# prod_all_subrun.sh on the cluster, one Slurm task per subrun.  Run by hand, mind
# that out/ here links to /data/cupsoft, which this machine cannot write to:
# give OUT=... instead.

set -u

DET=${1:-}
RUN=${2:-}
SUB=${3:-}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT=${OUT:-$HERE/out}

case "$DET" in
  fadc) exec "$HERE/prod_wcmd_subrun.sh" "$SUB" -r "$RUN" -o "$OUT" ${CFG:+-c "$CFG"} ${4:+"$4"} ;;
  iadc) exec "$HERE/prod_psmd_subrun.sh" "$SUB" -r "$RUN" -o "$OUT" ${CFG:+-c "$CFG"} ;;
  *) echo "usage: $0 <fadc|iadc> <run> <subrun> [boards]"; exit 2 ;;
esac
