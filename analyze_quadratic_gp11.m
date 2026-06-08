// analyze_quadratic_gp11.m -- Per-survivor analysis for the quadratic (1,11) search.
// Reads survivors_quadratic_gp11.m (from search_quadratic_gp11.m -- no re-inversion)
// and for each survivor: builds the curve over K, applies the minimal quadratic twist,
// reports conductor exponents at odd primes, and prints the Frobenius polynomial at 11.
//
// Direct analogue of analyze_survivors.m for the (1,5) case.
//
// Run: magma analyze_quadratic_gp11.m
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, ConductorExponentAt, PrimesAbove;

load "survivors_quadratic_gp11.m";  // defines k<s2> and survivors_quadratic_gp11
TwistBound := 30;
ZT<T> := PolynomialRing(Integers());
F11 := GF(11); R11 := PolynomialRing(F11);
p11 := PrimesAbove(k, 11)[1];

// Clear coefficient denominators for a curve over a number field K.
function IntegralModelK(C)
    K := BaseRing(C); f, h := HyperellipticPolynomials(C);
    m := LCM([Integers()| Denominator(K!c) : c in Coefficients(f) cat Coefficients(h) cat [K!1]]);
    return HyperellipticCurve(Parent(f)!(m^2*f), Parent(f)!(m*h));
end function;

// Reduce the integral curve CI modulo a prime ideal P of K.
function ReduceModP(CI, P)
    K := BaseRing(CI); f, h := HyperellipticPolynomials(CI);
    try
        Fp, redP := ResidueClassField(P); OK := Integers(K); Rx<x> := PolynomialRing(Fp);
        fp := Rx![redP(OK!c) : c in Coefficients(f)];
        hp := Rx![redP(OK!c) : c in Coefficients(h)];
        if LeadingCoefficient(fp) eq 0 then return false, _; end if;
        return true, HyperellipticCurve(fp, hp);
    catch e; return false, _; end try;
end function;

printf "Analyzing %o survivor(s).\n", #survivors_quadratic_gp11;
for s in survivors_quadratic_gp11 do
    ptors := s[1]; QI := s[2]; ht := s[3];
    printf "\n=== Ptors = %o   (height %o) ===\n", ptors, ht;

    ok := true;
    try C := Genus2CurveFromIgusa(QI, k);
    catch e; printf "  curve build failed: %o\n", e`Object; ok := false; end try;
    if not ok then continue; end if;
    printf "  C = %o\n", C;

    Cmin := C; dtwist := k!1;
    try
        Cmin, dtwist := MinimalTwist(C, k : Bound := TwistBound);
        good11 := &and[ConductorExponentAt(Cmin, p) eq 0 : p in PrimesAbove(k, 11)];
        printf "  minimal twist d=%o;  good at 11: %o\n", dtwist, good11;
    catch e; printf "  MinimalTwist: %o\n", e`Object; end try;

    badodd := [];
    for ell in PrimesUpTo(TwistBound) do
        if ell eq 2 then continue; end if;
        for p in PrimesAbove(k, ell) do
            e := ConductorExponentAt(Cmin, p);
            if e gt 0 then Append(~badodd, <Norm(p), e>); end if;
        end for;
    end for;
    printf "  bad odd primes (norm:exp, up to %o): %o\n", TwistBound, badodd;

    try
        ok11, Cp11 := ReduceModP(IntegralModelK(Cmin), p11);
        error if not ok11, "bad/non-integral reduction";
        fp := ZT ! LPolynomial(Cp11);
        printf "  Frobenius at p|11 (norm %o): %o = %o  (mod-11 roots: %o)\n",
            Norm(p11), fp, Factorization(fp), [r[1] : r in Roots(R11!fp)];
    catch e; printf "  Frobenius at 11: %o\n", e`Object; end try;
end for;
