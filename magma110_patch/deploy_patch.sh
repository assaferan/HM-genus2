#!/usr/bin/env bash
# Deploy the Magma #110 workaround as a PRIVATE, patched Magma copy.
#
# Touches no system files and redistributes no Magma source: it copies an
# existing Magma install you already have a licence for, then applies our
# ~25-line diff (definite.m.patch) to that copy.  Run the copy's binary for
# the affected computations; use stock Magma for everything else.
#
# Usage:  ./deploy_patch.sh <magma-install-dir> <dest-dir>
#   e.g.  ./deploy_patch.sh /opt/magma/magma-2.29-9 /scratch/$USER/magma-patched
#
# See README.md for the bug (Magma #110), the root cause, and the validation.
set -euo pipefail

SRC="${1:?usage: deploy_patch.sh <magma-install-dir> <dest-dir>}"
DEST="${2:?usage: deploy_patch.sh <magma-install-dir> <dest-dir>}"
HERE="$(cd "$(dirname "$0")" && pwd)"

[ -x "$SRC/magma" ] || { echo "error: no magma binary at $SRC/magma"; exit 1; }
[ -e "$DEST" ]      && { echo "error: dest $DEST already exists; remove it first"; exit 1; }

echo "1/3  copying $SRC -> $DEST (a few GB) ..."
cp -a "$SRC" "$DEST"

DEF="$DEST/package/Geometry/ModFrmHil/definite.m"
[ -f "$DEF" ] || { echo "error: definite.m not found at $DEF"; exit 1; }

echo "2/3  applying definite.m.patch ..."
patch "$DEF" < "$HERE/definite.m.patch"
rm -f "$DEST/package/Geometry/ModFrmHil/definite.sig"   # force recompile from patched .m

echo "3/3  done.  Patched Magma binary: $DEST/magma"
echo
echo "Sanity check (should print a marker-free dim + Hecke trace with no 'No solution' error):"
echo "  $DEST/magma <<'EOF'"
echo "  K<a>:=QuadraticField(2); OK:=Integers(K);"
echo "  M:=HilbertCuspForms(K, 79*OK, [2,2]); P:=Factorization(3*OK)[1][1];"
echo "  printf \"dim=%o tr=%o\\\\n\", Dimension(M), Trace(HeckeOperator(M,P)); exit;"
echo "  EOF"
