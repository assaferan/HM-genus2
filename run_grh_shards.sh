#!/bin/bash
# Launch grh_kernel.m per matched curve (idx 1..36; idx 37-39 blocked by Magma #110).
# Tiered BOUND: 3000 for dim <=9k, 800 for the 9 large-dim curves (still >> the FS ~200 bound).
cd /scratch/home/assaferan/GitHub/HM-genus2 || exit 1
MAGMA=${MAGMA:-/usr/local/bin/magma}
BIG="22 23 24 25 26 27 34 35 36"   # dim >9k: 243049, 312769x2, 328329, 478593x2, 377233, 472993x2
mkdir -p grhshards
echo "launching GRH certificate shards at $(date)"
for i in $(seq 1 36); do
  B=3000; for b in $BIG; do [ "$i" = "$b" ] && B=800; done
  nohup "$MAGMA" idx:="$i" BOUND:="$B" grh_kernel.m < /dev/null > "grhshards/$i.log" 2>&1 &
  echo "  shard $i (BOUND=$B) -> pid $!"
done
echo "all 36 launched; watch: grep -h 'GRH-KERNEL CERT' grh_kernel_*.out"
