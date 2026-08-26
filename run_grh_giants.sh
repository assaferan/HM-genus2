#!/bin/bash
# GRH certificates for the three Q(sqrt3) giants (idx 37/38/39), which need the
# Magma #110 workaround (see magma110_patch/): grh_kernel.m calls HilbertCuspForms,
# so on STOCK Magma <= V2.29-9 these crash at definite.m:1060 exactly like the match did.
# (Fixed upstream in V2.29-10 -- on such a host, point PATCHED at the stock binary.)
#
# SEQUENTIAL, always. Each of these spaces is dim 47728-76606 and runs 400+ GB; two at once
# do not fit on a 1 TB host, and running two in parallel is what destroyed the 785473 matches
# on 2026-08-18.
#
# Cost warning: unlike the match (14 fingerprint primes), the certificate verifies EVERY good
# prime N(P) <= BOUND -- one HeckeOperator each -- on top of re-deriving the match level, so each
# giant repeats the ~30 h space build and then adds the trace loop. Measured reference points:
# dim 39418 ran ~340 s/prime; dim 76606 is appreciably worse. Budget DAYS per giant.
# BOUND=400 is already well above the GRH Faltings-Serre bound ~(log cond)^2 ~ 185 here, and the
# trace loop is linear in the prime count, so do not raise it "just to be safe".
#
#   bash run_grh_giants.sh
#   REPO=~/GitHub/HM-genus2 PATCHED=~/magma-patched/magma IDXS="37 38" bash run_grh_giants.sh
#
# Writes grh_kernel_<idx>.out; sentinel per curve: "GRH-KERNEL CERT <label>: MODULAR".
set -u
cd "${REPO:-/scratch/home/assaferan/GitHub/HM-genus2}" || exit 1
PATCHED=${PATCHED:-/scratch/home/assaferan/magma-patched/magma}
BOUND=${BOUND:-400}
IDXS=${IDXS:-"37 38 39"}
FLOOR_GB=${FLOOR_GB:-700}
POLL=${POLL:-300}
[ -x "$PATCHED" ] || { echo "patched magma not found at $PATCHED"; exit 1; }
mkdir -p grhshards
avail() { awk '/^MemAvailable:/{print int($2/1048576)}' /proc/meminfo; }

echo "GRH giants on $PATCHED, BOUND=$BOUND, sequential, floor ${FLOOR_GB}G -- $(date)"
for i in $IDXS; do
  # Don't redo a cert we already have -- these cost days.
  if grep -q "GRH-KERNEL CERT" "grh_kernel_$i.out" 2>/dev/null; then
    echo "[$(date +%H:%M)] idx $i already certified, skipping: $(grep -h 'GRH-KERNEL CERT' grh_kernel_$i.out | tail -1)"
    continue
  fi
  # The cert re-derives the match, but a confirmed match is still the sane precondition.
  if ! grep -q "MATCH" "kernel_$i.out" 2>/dev/null; then
    echo "[$(date +%H:%M)] *** idx $i has no confirmed MATCH in kernel_$i.out -- skipping ***"
    continue
  fi

  while [ "$(avail)" -lt "$FLOOR_GB" ]; do
    echo "[$(date +%H:%M)] avail $(avail)G < ${FLOOR_GB}G -- holding before idx $i"
    sleep "$POLL"
  done

  echo "=== idx $i : starting $(date) (avail $(avail)G) ==="
  "$PATCHED" -b idx:="$i" BOUND:="$BOUND" grh_kernel.m < /dev/null > "grhshards/grh_giant_$i.log" 2>&1
  rc=$?

  if grep -q "GRH-KERNEL CERT" "grh_kernel_$i.out" 2>/dev/null; then
    echo "=== idx $i : done $(date) (rc=$rc) ==="
    grep -h "GRH-KERNEL CERT" "grh_kernel_$i.out" | tail -1
  else
    # Same guard as run_giants_sequential.sh: a vanished process is not a finished one.
    echo "=== idx $i : ENDED WITHOUT A CERT (rc=$rc, avail $(avail)G) ==="
    # Don't guess at the cause: rc and free RAM together identify it. rc=137 (SIGKILL) with low
    # avail is an OOM; rc=143 (SIGTERM) with plenty free is something/someone killing the job --
    # on a shared CI host, plausibly a cleanup step that reaps processes by name.
    echo "    No sentinel in grh_kernel_$i.out. rc=$rc ($([ $rc -eq 137 ] && echo SIGKILL; [ $rc -eq 143 ] && echo SIGTERM; [ $rc -ne 137 ] && [ $rc -ne 143 ] && echo "exit $rc")), $(avail)G free at exit."
    echo "    STOPPING rather than feeding another multi-day run into the same conditions."
    tail -3 "grhshards/grh_giant_$i.log" 2>/dev/null | sed 's/^/    | /'
    exit 1
  fi
done
echo "all requested giants done -- $(date)"
grep -h "GRH-KERNEL CERT" grh_kernel_3[789].out 2>/dev/null
