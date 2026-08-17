#!/bin/bash
# GRH certificates for the three Q(sqrt3) giants (idx 37/38/39), which need the
# Magma #110 workaround (see magma110_patch/): grh_kernel.m calls HilbertCuspForms,
# so on STOCK Magma these crash at definite.m:1060 exactly like the match did.
#
# Run these on the PATCHED Magma, and SEQUENTIALLY: each space is dim 47k-55k and
# ~150-200GB RSS, and unlike the match (14 primes) the cert verifies every good prime
# N(P) <= BOUND -- one HeckeOperator each -- so it is long. BOUND=400 is already well
# above the GRH Faltings-Serre bound ~(log cond)^2 ~ 185 for these conductors.
#
# Prereqs: matches confirmed (kernel_37/38/39.out show survivor=1/control=0) and RAM free.
# Writes grh_kernel_<idx>.out; sentinel per curve: "GRH-KERNEL CERT <label>: MODULAR".
set -u
cd /scratch/home/assaferan/GitHub/HM-genus2 || exit 1
PATCHED=${PATCHED:-/scratch/home/assaferan/magma-patched/magma}
BOUND=${BOUND:-400}
[ -x "$PATCHED" ] || { echo "patched magma not found at $PATCHED"; exit 1; }
mkdir -p grhshards
echo "GRH giants on PATCHED magma ($PATCHED), BOUND=$BOUND, sequential — $(date)"
for i in 37 38 39; do
  echo "=== idx $i : starting $(date) ==="
  "$PATCHED" idx:="$i" BOUND:="$BOUND" grh_kernel.m < /dev/null > "grhshards/grh_giant_$i.log" 2>&1
  echo "=== idx $i : done $(date) ==="
  grep -h "GRH-KERNEL CERT" "grh_kernel_$i.out" 2>/dev/null || echo "  (no cert line — inspect grhshards/grh_giant_$i.log)"
done
echo "all giants done — $(date)"

# --- Also outstanding: relaunch the 3 GRH giants killed for RAM (idx 34/35/36). Those use
# --- STOCK magma (they do NOT crash; they were only stopped mid trace-verify). Run e.g.:
#   for i in 34 35 36; do nohup /usr/local/bin/magma idx:=$i BOUND:=800 grh_kernel.m \
#       < /dev/null > grhshards/$i.log 2>&1 & done
# --- but keep total giant-sized jobs <= 2 at once (2TB box, no swap).
