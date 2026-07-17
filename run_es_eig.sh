#!/bin/bash
# run_es_eig.sh MAXN NCHUNK -- parallel, interleaved, checkpointed eigenvalue precompute.
# Merges the base file + all partials into es_eig/<label>.txt (dedup by ideal label).
MAXN=${1:-15000}
NCHUNK=${2:-8}
cd /scratch/home/assaferan/GitHub/HM-genus2
LABEL="2.2.8.1-5000.1-j"
mkdir -p es_eig
rm -f "es_eig/${LABEL}.txt.fast" es_eig/${LABEL}.part.* es_eig/chunk_*.log
echo "launching $NCHUNK interleaved chunks to norm $MAXN"
for i in $(seq 0 $((NCHUNK-1))); do
  magma MAXN:=$MAXN NCHUNK:=$NCHUNK IDX:=$i es_eig_chunk.m > "es_eig/chunk_$i.log" 2>&1 &
done
wait
echo "chunks done; merging base + partials (dedup)"
# base + partials, keep first occurrence of each ideal-label key
cat "es_eig/${LABEL}.txt" es_eig/${LABEL}.part.* 2>/dev/null | awk -F: '!seen[$1]++' > "es_eig/${LABEL}.merged"
mv "es_eig/${LABEL}.merged" "es_eig/${LABEL}.txt"
echo "merged lines: $(wc -l < es_eig/${LABEL}.txt)"
grep -h "chunk " es_eig/chunk_*.log 2>/dev/null || true
