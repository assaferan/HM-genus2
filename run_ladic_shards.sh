#!/bin/bash
# Launch ladic_degree.m per matched curve: certify [E:Q] > 1 from the mod-l generalized
# eigenspace, i.e. "sigma does not arise from an elliptic curve" (results §4f), WITHOUT
# isolating the eigenform.
#
# Tier 1 (default IDXS): the 27 matched forms of dim <= 8345 -- small enough to run
# alongside the §4d giant match jobs. Tier 2 = idx 22-27 (dim 11371-22719): run only
# once the giants have freed their RAM (IDXS="22 23 24 25 26 27").
#
# Unlike run_grh_shards.sh (which launched all 36 at once and contributed to an OOM),
# this caps concurrency and refuses to start a shard when free RAM is low.
#
#   bash run_ladic_shards.sh                       # tier 1, 6 at a time, 150G floor
#   IDXS="22 23 24" MAXJOBS=2 bash run_ladic_shards.sh
#
# Collect:  grep -h LADIC-DEG ladicshards/*.log
set -u
cd /scratch/home/assaferan/GitHub/HM-genus2 || exit 1
MAGMA=${MAGMA:-/usr/local/bin/magma}
IDXS=${IDXS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 28 29 30 31 32 33"}
MAXJOBS=${MAXJOBS:-6}
FLOOR_GB=${FLOOR_GB:-150}
mkdir -p ladicshards

avail() { free -g | awk '/^Mem:/{print $7}'; }

echo "launching ladic degree shards at $(date)  (MAXJOBS=$MAXJOBS, FLOOR=${FLOOR_GB}G)"
for i in $IDXS; do
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do sleep 30; done
  while [ "$(avail)" -lt "$FLOOR_GB" ]; do
    echo "  [$(date +%H:%M)] avail $(avail)G < ${FLOOR_GB}G -- holding before idx $i"
    sleep 300
  done
  nohup "$MAGMA" idx:="$i" ladic_degree.m < /dev/null > "ladicshards/$i.log" 2>&1 &
  echo "  shard $i -> pid $! (avail $(avail)G)"
  sleep 5
done
wait
echo "all shards done at $(date)"
grep -h LADIC-DEG ladicshards/*.log | sort
