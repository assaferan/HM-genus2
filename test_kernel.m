// test_kernel.m -- gate test for the kernel-intersection matcher (kernel_torsion.m).
// Reproduces the known match 14303.1 (Q(sqrt2), l=11, level norm 14303, e=0, dim 595)
// and asserts the core invariant: the mod-l Hecke eigen-kernel survivor is exactly
// 1-dimensional under the true fingerprint and collapses to 0 under the t->t+1 control.
// A regression in the fingerprint, HeckeOperator handling, or the kernel intersection
// breaks the assert (no sentinel) -> gate FAIL.  Self-contained (no CHIMP), ~20-30s.
//
// Success sentinel (for run_tests.sh):  KERNEL TEST: PASS
SetColumns(0);
load "torsion_data.m";

// pick 14303.1 out of the dataset (do not hard-code coeffs)
cand := [r : r in torsion_data | r[3] eq "14303.1"];
assert #cand eq 1;
row := cand[1];
d := row[1]; l := row[2]; fc := row[4]; hc := row[5]; cc := row[7];
assert d eq 2 and l eq 11;

K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K); Fl := GF(l);
f := R![K| c[1]+c[2]*a : c in fc]; hh := R![K| c[1]+c[2]*a : c in hc];
C := HyperellipticCurve(f, hh);
base := ideal<OK | cc[1]+cc[2]*a>;               // prime-to-2 conductor (given for Q(sqrt2))
assert Norm(base) eq 14303;
N := base;                                        // e = 0 for 14303

// fingerprint t_P = -#C(F_P) mod l at good primes P coprime to base, 2, l
fp := [];
for pp in PrimesUpTo(60) do
    if #fp ge 10 then break; end if;
    if pp eq 2 or pp eq l then continue; end if;
    for tup in Factorization(pp*OK) do
        P := tup[1]; if Norm(base+P) ne 1 then continue; end if;
        k,red := ResidueClassField(P); Rk := PolynomialRing(k);
        Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(f)], Rk![red(c):c in Coefficients(hh)]);
        Append(~fp, <P, Fl!(-#Ck)>);
        if #fp ge 10 then break; end if;
    end for;
end for;
assert #fp ge 6;

M := HilbertCuspForms(K, N, [2,2]); dm := Dimension(M);
printf "space dim = %o (expected 595)\n", dm;
assert dm eq 595;

survivor := function(tw)
    V := VectorSpace(Fl, dm); cur := V; I := IdentityMatrix(Fl, dm);
    for pr in fp do
        Tp := ChangeRing(Matrix(HeckeOperator(M, pr[1])), Fl);
        cur := cur meet Kernel(Tp - (pr[2]+tw)*I);
        if Dimension(cur) eq 0 then break; end if;
    end for;
    return Dimension(cur);
end function;

surv := survivor(0);
ctrl := survivor(1);
printf "survivor = %o (expected 1),  control = %o (expected 0)\n", surv, ctrl;
assert surv eq 1;      // true fingerprint isolates exactly the matching eigenform
assert ctrl eq 0;      // shifted fingerprint discriminates: no false positive

printf "KERNEL TEST: PASS (14303.1 survivor=1 control=0)\n";
exit;
