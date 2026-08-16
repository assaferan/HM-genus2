#!/bin/bash
# Launch one detached Magma process per curve (39 total) for the sharded HMF sweep.
# Each shard writes sweep_<idx>.out + hecke_cutters_<idx>.m and logs to shards/<idx>.log.
cd /scratch/home/assaferan/GitHub/HM-genus2 || exit 1
MAGMA=${MAGMA:-/usr/local/bin/magma}
N=$(grep -c '^Append(~torsion_data' torsion_data.m)
echo "launching $N shards at $(date)"
mkdir -p shards
for i in $(seq 1 "$N"); do
  nohup "$MAGMA" idx:="$i" sweep_shard.m < /dev/null > "shards/$i.log" 2>&1 &
  echo "  shard $i -> pid $!"
done
echo "all $N shards launched; watch with: grep -h . sweep_*.out"
