#!/bin/bash
# es_setup.sh -- reproduce the EichlerShimuraHMF setup for the B_f period computation.
# Clones ../EichlerShimuraHMF (sibling) and registers form 2.2.8.1-5000.1-j in its Labels.m.
# Requires ../CHIMP to already be checked out (same as the rest of this project).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
GH="$(dirname "$HERE")"          # parent dir holding CHIMP / EichlerShimuraHMF as siblings
cd "$GH"

if [ ! -d CHIMP ]; then
  echo "WARNING: ../CHIMP not found. Check out CHIMP as a sibling of this repo first." >&2
fi

if [ ! -d EichlerShimuraHMF ]; then
  echo "cloning EichlerShimuraHMF..."
  git clone https://github.com/edgarcosta/EichlerShimuraHMF
fi

LAB="EichlerShimuraHMF/src/Labels.m"
if grep -q "2.2.8.1-5000.1-j" "$LAB"; then
  echo "form already registered in $LAB"
else
  echo "registering 2.2.8.1-5000.1-j in $LAB ..."
  python3 - "$LAB" <<'PY'
import sys
p = sys.argv[1]
lines = open(p).read().splitlines(keepends=True)
entry = '  <"2.2.8.1-5000.1-j", [[0,50]], 2, [<"7.1", -1*2>], x^2 + 2*x - 2>, // B_f (mod-5 idx-5): abelian surface, RM by Q(sqrt3)\n'
out, done = [], False
for ln in lines:
    out.append(ln)
    if (not done) and ('2.2.8.1-2738.1-f' in ln):
        out.append(entry); done = True
assert done, "anchor line 2.2.8.1-2738.1-f not found in Labels.m"
open(p, 'w').write(''.join(out))
print("  inserted.")
PY
fi

echo "done. Attach with:"
echo '  AttachSpec("../CHIMP/CHIMP.spec"); AttachSpec("../EichlerShimuraHMF/src/spec");'
echo "Then: magma B:=60 MAXN:=<norm bound in es_eig/PROGRESS.md> goal2_mod5_es_tau.m"
