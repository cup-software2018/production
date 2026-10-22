#!/bin/bash
# Run the production over a run on the Slurm cluster, one job per subrun,
# each job producing FADC1, FADC2 and the IADC of its subrun in turn.
#
#   ./submit.sh <run> [njobs] [from] [to]
#
# Defaults: 8 at a time, every subrun with raw data for any detector.
#
#   ./submit.sh <run> 8              every subrun of the run
#   ./submit.sh <run> 8 0 600        subruns 0..600 only
#   CFG_FADC=scan.yaml ./submit.sh <run> 8 0 9
#
# CFG_FADC and CFG_IADC replace prod_fadc.yaml and prod_iadc.yaml for this
# submission.  Each Slurm task runs prod_all.sh for one subrun, which runs
# prod.sh for each detector.  One subrun per job rather than a contiguous
# block each: a slow subrun then holds up one slot instead of a whole block,
# and a job that dies costs one subrun to redo.  The output is one file per
# subrun per stream, which is what a TChain wants anyway, so nothing is merged
# afterwards.
#
# The jobs run on qnode1-2 behind qserver (192.168.0.160), submitted over ssh.
# The nodes see only qserver's /home and /data, not this machine's /scratch,
# and cupsoft is uid 1000 here but 1005 there, so nothing here can write to
# /data/cupsoft.  The programs are built here (the nodes have no yaml-cpp
# headers) and copied over through qserver, with the scripts and the configs,
# into a fresh stage directory per submission; the output lands in $RDIR/out,
# which out/ here links to.  A resubmission never touches the programs under
# jobs still running from the last one.
#
# The script waits for the jobs and then summarises.  Ctrl-C stops the
# waiting only; the jobs keep running (squeue / scancel on qserver).

set -u

RUN=${1:-}
NJOBS=${2:-8}
FROM=${3:--1}
TO=${4:--1}
if ! [[ "$RUN" =~ ^[0-9]+$ && "$NJOBS" =~ ^[0-9]+$ && "$FROM" =~ ^-?[0-9]+$ && "$TO" =~ ^-?[0-9]+$ ]]; then
  echo "usage: $0 <run> [njobs] [from] [to]"; exit 2
fi

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CFG_FADC=$(readlink -f "${CFG_FADC:-$HERE/prod_fadc.yaml}")
CFG_IADC=$(readlink -f "${CFG_IADC:-$HERE/prod_iadc.yaml}")
RUNSTR=$(printf "%06d" "$RUN")
DATA=/data/amore2/data/RAW/$RUNSTR

QHOST=${QHOST:-cupsoft@192.168.0.160}
RDIR_DEFAULT=/data/cupsoft/WCMD/production
RDIR=${RDIR:-$RDIR_DEFAULT}
OUT=$RDIR/out
STAGE=$RDIR/stage/${RUNSTR}_$(date +%Y%m%d_%H%M%S)
MAXARRAY=${MAXARRAY:-1000}  # the cluster's MaxArraySize is 1001
q() { ssh -o BatchMode=yes "$QHOST" "$@"; }

# thisroot.sh appends to these, so they must exist under `set -u`.
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH:-}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:-}

source /home/cupsoft/Products/ROOT/V6.32.20_3.11.14/bin/thisroot.sh

[ -d "$DATA" ] || { echo "no raw directory $DATA"; exit 1; }
for c in "$CFG_FADC" "$CFG_IADC"; do [ -f "$c" ] || { echo "no config $c"; exit 1; }; done

# Build here; the nodes have no yaml-cpp headers to compile against.
{ cmake -S "$HERE" -B "$HERE/build" && cmake --build "$HERE/build"; } \
  >"$HERE/build.log" 2>&1 || { echo "build failed, see $HERE/build.log"; exit 1; }

# The subruns with raw data for at least one stream, one per line, and the
# log prod.sh will write for each stream that has it.
list=$(mktemp)
logs=$(mktemp)
trap 'rm -f "$list" "$logs"' EXIT
nstream=0
for prefix in FADC1 FADC2 IADC; do
  if [ "$FROM" -ge 0 ]; then
    subs=$(seq "$FROM" "$([ "$TO" -ge "$FROM" ] && echo "$TO" || echo "$FROM")")
  else
    subs=$(for f in "$DATA"/${prefix}_${RUNSTR}.root.*; do
             [ -e "$f" ] && echo $((10#${f##*.}))
           done)
  fi
  n=0
  for sub in $subs; do
    sub5=$(printf "%05d" "$sub")
    [ -e "$DATA/${prefix}_${RUNSTR}.root.$sub5" ] || continue
    echo "$sub" >>"$list"
    echo "$OUT/log/${prefix}_$sub5.log" >>"$logs"
    n=$((n + 1))
  done
  echo "$prefix: $n subruns"
  [ "$n" -gt 0 ] && nstream=$((nstream + 1))
done
[ "$nstream" -gt 0 ] || { echo "no input in $DATA"; exit 1; }
sort -n -u -o "$list" "$list"
nsub=$(wc -l <"$list")

q "mkdir -p '$STAGE/slurm' '$OUT/log'" || {
  echo "cannot reach $QHOST"; exit 1; }
# prod.sh finds produce_<det> and prod_<det>.yaml beside itself.
rsync -a "$HERE/build/produce_fadc" "$HERE/build/produce_iadc" "$HERE/prod.sh" \
  "$HERE/prod_all.sh" "$QHOST:$STAGE/" &&
rsync -a "$CFG_FADC" "$QHOST:$STAGE/prod_fadc.yaml" &&
rsync -a "$CFG_IADC" "$QHOST:$STAGE/prod_iadc.yaml" &&
rsync -a "$list" "$QHOST:$STAGE/subruns.txt" || {
  echo "staging to $QHOST:$STAGE failed"; exit 1; }

# out/ follows the default RDIR only, so a test with RDIR=... does not move
# the link the analysis code reads through.
if [ "$RDIR" != "$RDIR_DEFAULT" ]; then
  :
elif [ -L "$HERE/out" ] || [ ! -e "$HERE/out" ]; then
  ln -sfn "$OUT" "$HERE/out"
else
  echo "note: $HERE/out is a real directory, not linked to $OUT"
fi

# Task i runs prod_all.sh on line OFFSET+i of subruns.txt.  $sub and the task id
# must be expanded on the node, not by the shell on qserver, so the command is
# quoted once more (%q) for the ssh hop.
task='sub=$(sed -n "$((OFFSET + SLURM_ARRAY_TASK_ID))p" '"$STAGE/subruns.txt"') && '
task+="OUT=$OUT exec $STAGE/prod_all.sh $RUN "'"$sub"'

# Arrays of at most MAXARRAY tasks, each starting after the one before it, so
# no more than NJOBS run at once however long the list is.
ids=""
dep=""
for ((off = 0; off < nsub; off += MAXARRAY)); do
  n=$(( nsub - off < MAXARRAY ? nsub - off : MAXARRAY ))
  id=$(q "sbatch --parsable --job-name=prod_$RUN --array=1-$n%$NJOBS \
            --cpus-per-task=1 --mem=2G $dep --export=ALL,OFFSET=$off \
            --output='$STAGE/slurm/%A_%a.out' --wrap $(printf %q "$task")") || {
    echo "sbatch failed"; exit 1; }
  ids="$ids${ids:+,}$id"
  dep="--dependency=afterany:$id"
done

echo "run $RUN, configs $CFG_FADC $CFG_IADC, $NJOBS jobs at a time"
echo "$nsub subrun jobs as job $ids on $QHOST"
echo "stage  $STAGE"
echo "output $OUT"
echo "started $(date '+%F %T')"
start=$(date +%s)

while :; do
  left=$(q "squeue -h -r -j $ids -o %T" 2>/dev/null | sort | uniq -c | tr -s ' \n' ' ')
  [ -n "$left" ] || break
  printf '\r  %s  %-60s' "$(date '+%T')" "$left"
  sleep 30
done
echo

el=$(( $(date +%s) - start ))
echo "--- done $(date '+%F %T'); $nsub subrun jobs in $((el/60))m $((el%60))s ---"
# Checked on qserver: this side's NFS view of /data can lag behind the nodes.
echo "output: $(q "ls '$OUT'/prod_*_${RUNSTR}.*.root 2>/dev/null | wc -l") of $(wc -l <"$logs") files in $OUT"
failed=$(q 'while read -r l; do grep -q "^exit 0 " "$l" 2>/dev/null || echo "$l"; done' <"$logs")
bad=$(q 'xargs grep -lE "ERROR|WARNING|Error in" 2>/dev/null' <"$logs")
if [ -n "$failed" ]; then
  echo "failed or never ran:"
  echo "$failed" | sed 's|^|  |'
fi
if [ -n "$bad" ]; then
  echo "ERROR, WARNING or ROOT error in:"
  echo "$bad" | sed 's|^|  |'
fi
