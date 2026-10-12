#!/bin/bash
# Produce one subrun of every detector, in turn: FADC1, FADC2, then the IADC.
#
#   ./prod_all.sh <run> <subrun>
#
# Each detector goes through prod.sh, with its default settings
# (prod_fadc.yaml, prod_iadc.yaml beside it) and the same RAW, OUT and BIN
# from the environment; see prod.sh.  A detector with no raw file for this
# subrun is skipped.  Exits non-zero if any producer fails or no detector has
# input.  submit.sh runs this on the cluster, one Slurm task per subrun.

set -u

RUN=${1:-}
SUB=${2:-}
if ! [[ "$RUN" =~ ^[0-9]+$ && "$SUB" =~ ^[0-9]+$ ]]; then
  echo "usage: $0 <run> <subrun>"; exit 2
fi
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# prod.sh exits 3 when none of its streams has a raw file for the subrun.
status=0
ran=0
for det in fadc iadc; do
  CFG=$HERE/prod_$det.yaml "$HERE/prod.sh" "$det" "$RUN" "$SUB"
  rc=$?
  case $rc in
    0) ran=$((ran + 1)) ;;
    3) ;;
    *) ran=$((ran + 1)); status=1 ;;
  esac
done
[ "$ran" -gt 0 ] || { echo "no raw input for run $RUN subrun $SUB"; exit 3; }
exit $status
