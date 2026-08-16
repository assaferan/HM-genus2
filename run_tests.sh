#!/usr/bin/env bash
# run_tests.sh -- the project's test gate (for CI and local use).
#
# Magma NEVER sets a non-zero exit code -- not even on a runtime error -- so this
# script cannot rely on Magma's exit status. Instead each Magma gate must print a
# unique SUCCESS SENTINEL on the happy path; a gate is considered PASS iff its
# combined stdout+stderr contains that sentinel AND contains no Magma error marker.
#
# Usage:   ./run_tests.sh
# Env:     RUN_CHIMP_TESTS=0   skip the CHIMP-dependent full-pipeline test even if ../CHIMP exists
#          MAGMA=/path/to/magma   override the Magma binary
set -u
cd "$(dirname "$0")"

MAGMA="${MAGMA:-magma}"
LOGDIR="$(mktemp -d)"
FAILURES=0
RAN=0

# Magma error markers that must never appear in a passing run.
ERR_RE='Runtime error|error, line|>> *$|Segmentation fault|Internal error|Assertion|Aborted'

pass() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; FAILURES=$((FAILURES+1)); }
skip() { printf '  \033[33mSKIP\033[0m  %s\n' "$1"; }

# run_magma <label> <sentinel> <magma-args...>
run_magma() {
  local label="$1" sentinel="$2"; shift 2
  local log="$LOGDIR/$(echo "$label" | tr -c 'A-Za-z0-9' '_').log"
  RAN=$((RAN+1))
  # -b: batch mode (no banner/prompt); </dev/null: never block waiting on stdin.
  "$MAGMA" -b "$@" </dev/null >"$log" 2>&1
  if grep -qiE "$ERR_RE" "$log"; then
    fail "$label (Magma error)"; sed 's/^/      | /' "$log" | grep -iE "$ERR_RE" | head -3; return
  fi
  if ! grep -qF "$sentinel" "$log"; then
    fail "$label (missing sentinel: \"$sentinel\")"; tail -3 "$log" | sed 's/^/      | /'; return
  fi
  pass "$label"
}

echo "== Environment =="
if ! command -v "$MAGMA" >/dev/null 2>&1; then
  echo "  Magma not found (MAGMA=$MAGMA). Cannot run the Magma gate." >&2
  exit 2
fi
echo "  magma: $(command -v "$MAGMA")"
echo "  logs:  $LOGDIR"

echo
echo "== Python (compile only) =="
if command -v python3 >/dev/null 2>&1; then
  py_fail=0
  for f in *.py; do
    [ -e "$f" ] || continue
    RAN=$((RAN+1))
    if python3 -m py_compile "$f" 2>"$LOGDIR/py.log"; then pass "py_compile $f"
    else fail "py_compile $f"; sed 's/^/      | /' "$LOGDIR/py.log"; py_fail=1; fi
  done
  [ "$py_fail" = 0 ] || true
else
  skip "python3 not found -- skipping py_compile"
fi

echo
echo "== Magma: self-contained gate (no CHIMP) =="
# Dataset structure + conductor-norm validation over all 39 curves.
run_magma "validate.m (dataset + 1+chi+sigma structure)" "VALIDATE: ALL PASS" validate.m
# Kernel-intersection matcher: reproduce the known 14303.1 match (survivor=1, control=0).
run_magma "test_kernel.m (kernel matcher: 14303.1 survivor/control)" "KERNEL TEST: PASS" test_kernel.m
# GRH certificate from the survivor eigenvector: reproduce 14303.1 MODULAR (small bound).
run_magma "grh_kernel.m (14303.1 GRH cert, BOUND=60)" "GRH-KERNEL CERT 14303.1: MODULAR" idx:=3 BOUND:=60 grh_kernel.m

echo
echo "== Magma: full-pipeline gate (needs CHIMP) =="
if [ "${RUN_CHIMP_TESTS:-1}" != "1" ]; then
  skip "RUN_CHIMP_TESTS=0 -- full-pipeline test disabled"
elif [ -f "../CHIMP/CHIMP.spec" ]; then
  # test.m needs Genus2Curve in scope; it prints "PASS: ..." and errors otherwise.
  run_magma "test.m (x in P^3 -> genus-2 curve, 5-isogeny)" "PASS: Jac(C) admits a rational cyclic 5-isogeny" Genus2Curve.m test.m
else
  skip "../CHIMP/CHIMP.spec not found -- full-pipeline test needs CHIMP as a sibling checkout"
fi

echo
echo "== Summary =="
echo "  ran $RAN gate(s), $FAILURES failure(s)"
if [ "$FAILURES" -ne 0 ]; then
  echo "  RESULT: FAIL"; exit 1
fi
echo "  RESULT: PASS"; exit 0
