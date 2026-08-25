#!/bin/bash
# Run the §4d giant matches ONE AT A TIME, waiting for each to finish before starting the next.
#
# Why sequential: each giant peaks somewhere north of 450 GB (on lovelace they were at 424-452 GB
# while still *building* when they were killed). Two at once fit in lovelace's 2 TB only barely,
# and do NOT fit in hensel's 1 TB. Running two in parallel is what lost both runs on 2026-08-18.
#
# The important guard: a giant's process disappearing does NOT mean it succeeded -- it may have been
# OOM-killed, which is exactly what happened before. We only advance to the next giant if the
# previous one wrote its completion sentinel. If a giant vanishes without one, we STOP and say so,
# rather than feeding another 40 h run into whatever killed the last one.
#
#   bash run_giants_sequential.sh                       # idx 38 then 39, patched magma
#   IDXS="39" MAGMA=~/magma-patched/magma bash run_giants_sequential.sh
#
# Progress: giants_sequential.log; per-run kernel_<idx>.out and kernel_<idx>_run.log
set -u
cd "${REPO:-$HOME/GitHub/HM-genus2}" || exit 1
MAGMA=${MAGMA:-$HOME/magma-patched/magma}
IDXS=${IDXS:-"38 39"}
FLOOR_GB=${FLOOR_GB:-700}     # a giant needs ~500-600 G; refuse to start without real headroom
POLL=${POLL:-300}

[ -x "$MAGMA" ] || { echo "no magma at $MAGMA"; exit 1; }
avail() { awk '/^MemAvailable:/{print int($2/1048576)}' /proc/meminfo; }

echo "=== giants sequential: [$IDXS] with $MAGMA, floor ${FLOOR_GB}G -- $(date) ==="

for i in $IDXS; do
  # Already done from an earlier run? Skip rather than redo 40 h.
  if grep -q "KERNEL SHARD $i DONE" "kernel_$i.out" 2>/dev/null; then
    echo "[$(date +%H:%M)] idx $i already complete, skipping:"
    grep -h "MATCH\|survivor" "kernel_$i.out" | tail -2
    continue
  fi

  while [ "$(avail)" -lt "$FLOOR_GB" ]; do
    echo "[$(date +%H:%M)] avail $(avail)G < ${FLOOR_GB}G -- holding before idx $i"
    sleep "$POLL"
  done

  echo "[$(date +%H:%M)] starting idx $i (avail $(avail)G)"
  "$MAGMA" -b idx:="$i" kernel_torsion.m < /dev/null > "kernel_${i}_run.log" 2>&1
  rc=$?

  if grep -q "KERNEL SHARD $i DONE" "kernel_$i.out" 2>/dev/null; then
    echo "[$(date +%H:%M)] idx $i COMPLETED (rc=$rc):"
    grep -h "MATCH\|survivor" "kernel_$i.out" | tail -2
  else
    echo "[$(date +%H:%M)] *** idx $i ENDED WITHOUT COMPLETING (rc=$rc, avail $(avail)G) ***"
    echo "    No completion sentinel in kernel_$i.out -- most likely OOM-killed."
    echo "    STOPPING: not starting the remaining giants until this is understood."
    tail -3 "kernel_${i}_run.log" 2>/dev/null | sed 's/^/    | /'
    exit 1
  fi
done

echo "=== all requested giants done -- $(date) ==="
grep -h "MATCH" kernel_3[789].out 2>/dev/null
