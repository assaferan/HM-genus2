#!/bin/bash
# Verify the banked evidence in evidence/ against torsion_data.m.
#
# The evidence is the actual product of weeks of compute, and the runs that produced it are
# long gone -- so it needs to be checkable without re-running anything. This is pure text
# checking: no Magma, runs in well under a second.
#
# It exists because of two near-misses worth guarding against permanently:
#   * kernel_37/38/39.out on the original host are stale HEADERS from runs that were OOM-killed
#     before writing a result. They look plausible -- right name, right size-ish -- but contain
#     no MATCH. Banking one of those would silently misrepresent an unfinished run as a result.
#   * every curve has a conjugate (881.1 / 881.2, ...) with near-identical output, so a
#     transposed file would be easy to miss by eye.
# Both are caught by checking the LABEL inside each file against torsion_data.m by index.
#
# Sentinel on success: "EVIDENCE: ALL PASS (...)"
set -u
cd "$(dirname "$0")" || exit 1
EV=${EV:-evidence}
fail=0
note() { echo "  $*"; }
bad()  { echo "  FAIL: $*"; fail=1; }

[ -d "$EV" ] || { echo "EVIDENCE: no $EV/ directory"; exit 1; }

# Labels in Append order == idx order. (Read this way, not with mapfile: macOS ships bash 3.2,
# and this script has to run on both the laptop and the Linux CI runner.)
LABELS=()
while IFS= read -r l; do LABELS+=("$l"); done < <(grep -o '"[0-9]*\.[0-9]*"' torsion_data.m | tr -d '"')
n=${#LABELS[@]}
if [ "$n" -ne 39 ]; then bad "expected 39 labels in torsion_data.m, found $n"; fi

# ---- matches: all 39 must be present, complete, and correctly labelled ----
nmatch=0
for i in $(seq 1 "$n"); do
  lab="${LABELS[$((i-1))]}"
  f="$EV/matches/kernel_$i.out"
  if [ ! -f "$f" ]; then bad "missing match evidence for idx $i ($lab)"; continue; fi
  grep -q "KERNEL SHARD $i DONE" "$f" || { bad "idx $i ($lab): no completion sentinel (unfinished/killed run?)"; continue; }
  line=$(grep -E "MATCH \(e=" "$f" | tail -1)
  [ -n "$line" ]                  || { bad "idx $i ($lab): no MATCH line"; continue; }
  case "$line" in "$lab"*) ;; *)     bad "idx $i: label mismatch -- expected '$lab', file says '${line%% *}'"; continue;; esac
  echo "$line" | grep -q "survivor=1 control=0" || { bad "idx $i ($lab): not survivor=1/control=0"; continue; }
  nmatch=$((nmatch+1))
done

# ---- certificates: however many exist must be MODULAR, 0 disagreements, right label ----
ncert=0; missing=""
for i in $(seq 1 "$n"); do
  lab="${LABELS[$((i-1))]}"
  f="$EV/certificates/grh_kernel_$i.out"
  if [ ! -f "$f" ]; then missing="$missing $i:$lab"; continue; fi
  line=$(grep "GRH-KERNEL CERT" "$f" | tail -1)
  [ -n "$line" ]                                  || { bad "cert idx $i ($lab): file present but no CERT line"; continue; }
  echo "$line" | grep -q "GRH-KERNEL CERT $lab: MODULAR" || { bad "cert idx $i: expected '$lab: MODULAR', got '$line'"; continue; }
  echo "$line" | grep -q "disagree=0"             || { bad "cert idx $i ($lab): disagreements present"; continue; }
  echo "$line" | grep -q "irred"                  || { bad "cert idx $i ($lab): sigma not recorded irreducible"; continue; }
  ncert=$((ncert+1))
done

note "matches:      $nmatch/$n complete, control-validated, labels agree with torsion_data.m"
note "certificates: $ncert/$n MODULAR, 0 disagreements, sigma irreducible"
if [ -n "$missing" ]; then note "not yet certified:$missing"; fi

# Evidence present on disk but not committed is evidence that does not exist for anyone else.
# .gitignore carries *.out and *.log rules for LaTeX artifacts, so `git add evidence/` silently
# skipped all 104 files the first time and shipped an empty directory -- green locally, red in CI.
# Check it here so that failure surfaces before the push, not after.
if [ "${EV}" = "evidence" ] && command -v git >/dev/null 2>&1 && git rev-parse --git-dir >/dev/null 2>&1; then
  untracked=0
  # Every evidence subdirectory, not just the three checked above -- an untracked file in a
  # directory this loop forgets is exactly the blind spot that shipped an empty evidence/ once.
  for f in "$EV"/matches/*.out "$EV"/certificates/*.out "$EV"/ladic/*.log "$EV"/image/*.txt; do
    [ -e "$f" ] || continue
    git ls-files --error-unmatch "$f" >/dev/null 2>&1 || untracked=$((untracked+1))
  done
  if [ "$untracked" -ne 0 ]; then
    bad "$untracked evidence file(s) exist on disk but are NOT tracked by git (check .gitignore)"
  else
    note "all evidence files are tracked by git"
  fi
fi

if [ "$fail" -ne 0 ]; then echo "EVIDENCE: FAIL"; exit 1; fi
if [ "$nmatch" -ne "$n" ]; then echo "EVIDENCE: FAIL (only $nmatch/$n matches)"; exit 1; fi
echo "EVIDENCE: ALL PASS ($nmatch/$n matches, $ncert/$n certificates)"
