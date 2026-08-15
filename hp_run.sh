#!/bin/bash
# hp_run.sh -- background high-precision eigenvalue precompute for B_f's period matrix.
# Loops over increasing norm targets; after each, merges the cache and commits + pushes
# es_eig/<label>.txt (+ PROGRESS.md + all session scripts) to branch `es-highprec`, so
# another machine can `git fetch && git checkout es-highprec` and read the results.
#
# Launch DETACHED so it survives the session:  nohup bash hp_run.sh >/dev/null 2>&1 &
cd /scratch/home/assaferan/GitHub/HM-genus2
LABEL="2.2.8.1-5000.1-j"
BRANCH="es-highprec"
NCHUNK=8
LOCK=".hp_run.lock"

# single-instance lock
if [ -e "$LOCK" ] && kill -0 "$(cat "$LOCK")" 2>/dev/null; then
  echo "hp_run already running (pid $(cat "$LOCK"))"; exit 0
fi
echo $$ > "$LOCK"
trap 'rm -f "$LOCK"' EXIT

exec >> hp_run.log 2>&1
echo "=========================================================="
echo "=== hp_run started $(date) pid $$ ==="

TRAILER=$'\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_017pZCVUsUCZG9naw89Xar8S'

commit_push () {
  git add -A
  git commit -m "$1$TRAILER" || echo "(nothing to commit)"
  git push -u origin "$BRANCH" || echo "push failed $(date) (will retry next milestone)"
}

# switch to the es-highprec branch, keeping the working tree
git rev-parse --verify "$BRANCH" >/dev/null 2>&1 || git branch "$BRANCH"
git checkout "$BRANCH"

# initial snapshot
{
  echo "# B_f eigenvalue precompute progress"
  echo ""
  echo "- started $(date)"
  echo "- baseline: contiguous norm bound 10103 (period matrix tau computed to ~4 digits)"
} > es_eig/PROGRESS.md
commit_push "ES handoff: scripts, notes, HANDOFF_ES.md, period matrix (~4 digits), eigenvalues to norm 10103"

for N in 20000 30000 50000 75000 100000 150000; do
  echo "--- target norm $N  $(date) ---"
  bash run_es_eig.sh "$N" "$NCHUNK" || echo "run_es_eig.sh returned nonzero for $N"
  BOUND=$(magma -b es_eig_bound.m 2>/dev/null | grep BOUND | cut -d: -f2)
  LINES=$(wc -l < "es_eig/$LABEL.txt")
  {
    echo "- $(date): target norm $N done; contiguous bound = ${BOUND:-?}; cache lines = $LINES"
    echo "    -> recompute period matrix at higher precision:  magma B:=60 MAXN:=${BOUND:-$N} goal2_mod5_es_tau.m"
  } >> es_eig/PROGRESS.md
  commit_push "ES eigenvalues: reached target norm $N (contiguous bound ${BOUND:-?}, $LINES lines)"
done

echo "=== hp_run finished $(date) ==="
echo "- $(date): high-precision run COMPLETE" >> es_eig/PROGRESS.md
commit_push "ES eigenvalues: high-precision run complete"
