#!/bin/bash
# Launch one Magma process per curve (39) for the kernel-intersection HMF matcher.
# Each shard writes kernel_<idx>.out and logs to kshards/<idx>.log.
cd /scratch/home/assaferan/GitHub/HM-genus2 || exit 1
MAGMA=${MAGMA:-/usr/local/bin/magma}
N=$(grep -c '^Append(~torsion_data' torsion_data.m)
echo "launching $N kernel shards at $(date)"
mkdir -p kshards
for i in $(seq 1 "$N"); do
  nohup "$MAGMA" idx:="$i" kernel_torsion.m < /dev/null > "kshards/$i.log" 2>&1 &
  echo "  kernel shard $i -> pid $!"
done
echo "all $N launched; watch: grep -h MATCH kernel_*.out"
