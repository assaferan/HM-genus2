#!/bin/bash
# Wait for the two §4d giant match jobs (idx 38/39) to exit AND for RAM to come back,
# then relaunch the three GRH certificate shards that were killed mid-trace-verify
# (idx 34/35/36 = 377233.2, 472993.1, 472993.2). Those run on STOCK magma -- they never
# crashed, they were only stopped to protect the giants' RAM.
#
# Guards, because this starts ~200-290GB jobs unattended:
#   - never starts while either giant PID is alive
#   - never exceeds MAXJOBS concurrent shards (default 2; the 2TB box has no swap)
#   - never starts a shard unless FLOOR_GB of RAM is free right now
#   - gives up waiting after TIMEOUT_H hours rather than lurking forever
#
#   bash watch_and_relaunch_grh.sh                 # defaults
#   GIANT_PIDS="3711411 3711412" FLOOR_GB=350 bash watch_and_relaunch_grh.sh
#
# Progress -> grh_watcher.log; certs -> grh_kernel_<idx>.out ("GRH-KERNEL CERT ...: MODULAR").
set -u
cd /scratch/home/assaferan/GitHub/HM-genus2 || exit 1
MAGMA=${MAGMA:-/usr/local/bin/magma}
GIANT_PIDS=${GIANT_PIDS:-"3711411 3711412"}
IDXS=${IDXS:-"34 35 36"}
BOUND=${BOUND:-800}
MAXJOBS=${MAXJOBS:-2}
FLOOR_GB=${FLOOR_GB:-350}
TIMEOUT_H=${TIMEOUT_H:-72}
mkdir -p grhshards

avail() { free -g | awk '/^Mem:/{print $7}'; }
alive() { for p in $GIANT_PIDS; do kill -0 "$p" 2>/dev/null && return 0; done; return 1; }

echo "watcher started $(date): waiting on giants [$GIANT_PIDS], then idx [$IDXS] at BOUND=$BOUND"
deadline=$(( $(date +%s) + TIMEOUT_H*3600 ))

while alive; do
  [ "$(date +%s)" -ge "$deadline" ] && { echo "TIMEOUT after ${TIMEOUT_H}h with giants still alive -- NOT launching."; exit 1; }
  sleep 600
done
echo "giants exited at $(date); match results:"
grep -h "MATCH\|survivor" kernel_38.out kernel_39.out 2>/dev/null | tail -4
grep -h "MATCH\|survivor" /scratch/home/assaferan/patch_validate/kernel_3[89].out 2>/dev/null | tail -4

for i in $IDXS; do
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do sleep 60; done
  while [ "$(avail)" -lt "$FLOOR_GB" ]; do
    echo "  [$(date +%H:%M)] avail $(avail)G < ${FLOOR_GB}G -- holding before idx $i"
    sleep 300
  done
  nohup "$MAGMA" idx:="$i" BOUND:="$BOUND" grh_kernel.m < /dev/null > "grhshards/$i.log" 2>&1 &
  echo "  idx $i -> pid $! at $(date +%H:%M) (avail $(avail)G)"
  sleep 30
done
wait
echo "all relaunched shards finished $(date)"
grep -h "GRH-KERNEL CERT" grh_kernel_3[456].out 2>/dev/null
