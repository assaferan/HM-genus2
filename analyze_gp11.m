// analyze_gp11.m -- Heavy per-survivor analysis for the (1,11) pipeline.
// Reads survivors_gp11.m (written by search_gp11.m -- NO period re-computation),
// and for each survivor: rebuilds K and the Igusa invariants, constructs the curve,
// takes the minimal quadratic twist (good reduction at primes of potential-good
// reduction at 11), reports conductor exponents at odd primes, and prints the
// factored Frobenius polynomial at (a prime above) 11.
//
// Run: magma analyze_gp11.m
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, ConductorExponentAt, PrimesAbove;

load "survivors_gp11.m";   // defines survivors_gp11
TwistBound := 30;

Qx<t_> := PolynomialRing(Rationals());
ZT<T>  := PolynomialRing(Integers());
F11 := GF(11); R11 := PolynomialRing(F11);

// Clear coefficient denominators for a curve over Q or a number field K.
function IntegralModelK(C)
    K := BaseRing(C);
    f, h := HyperellipticPolynomials(C);
    if Type(K) eq FldRat then
        m := LCM([Integers()| Denominator(c) : c in Coefficients(f) cat Coefficients(h) cat [1]]);
        return HyperellipticCurve(Parent(f)!(m^2*f), Parent(f)!(m*h));
    end if;
    m := LCM([Integers()| Denominator(K!c) : c in Coefficients(f) cat Coefficients(h) cat [K!1]]);
    return HyperellipticCurve(Parent(f)!(m^2*f), Parent(f)!(m*h));
end function;

// Reduce the integral curve CI modulo prime P (integer for K=Q, prime ideal otherwise).
function ReduceModP(CI, P)
    K := BaseRing(CI); f, h := HyperellipticPolynomials(CI);
    try
        if Type(K) eq FldRat then
            Fp := GF(P); Rx<x> := PolynomialRing(Fp);
            fp := Rx![Fp!Integers()!c : c in Coefficients(f)];
            hp := Rx![Fp!Integers()!c : c in Coefficients(h)];
        else
            Fp, redP := ResidueClassField(P);
            OK := Integers(K); Rx<x> := PolynomialRing(Fp);
            fp := Rx![redP(OK!c) : c in Coefficients(f)];
            hp := Rx![redP(OK!c) : c in Coefficients(h)];
        end if;
        if LeadingCoefficient(fp) eq 0 then return false, _; end if;
        return true, HyperellipticCurve(fp, hp);
    catch e; return false, _; end try;
end function;

// Frobenius (L-)polynomial of C at P as a ZT-polynomial.
function FrobPolyAt(C, P)
    ok, Cp := ReduceModP(IntegralModelK(C), P);
    error if not ok, "bad/non-integral reduction at P";
    return ZT ! LPolynomial(Cp);
end function;

// True iff C has actual good reduction at all primes of K above 11.
function GoodReductionAt11(C)
    return &and[ConductorExponentAt(C, p) eq 0 : p in PrimesAbove(BaseRing(C), 11)];
end function;

printf "Analyzing %o survivor(s) from survivors_gp11.m.\n\n", #survivors_gp11;

for s in survivors_gp11 do
    a := s[1]; h_coeffs := s[2]; QI_seqs := s[3]; ht := s[4];

    // Rebuild K from h's coefficient list (degree 0 upward) and QI from coord-sequences.
    h := Qx ! h_coeffs;    // Qx ! [c0,...,cd]  ->  c0 + c1*t_ + ... + cd*t_^d
    d := Degree(h);
    if d eq 1 then
        K := Rationals();
        QI := [Rationals() ! qi[1] : qi in QI_seqs];
    else
        K<aa> := NumberField(h);
        QI := [K ! qi : qi in QI_seqs];
    end if;
    Kstr := Type(K) eq FldRat select "Q" else Sprint(DefiningPolynomial(K));
    printf "=== a=%o  ht=%o  K=%o ===\n", a, ht, Kstr;

    ok := true;
    try C := Genus2CurveFromIgusa(QI, K);
    catch e; printf "  curve build failed: %o\n", e`Object; ok := false; end try;
    if not ok then continue; end if;
    printf "  C = %o\n", C;

    // Minimal quadratic twist (targets good reduction at all odd primes of bad reduction).
    Cmin := C; dtwist := Type(K) eq FldRat select 1 else K!1;
    try
        Cmin, dtwist := MinimalTwist(C, K : Bound := TwistBound);
        printf "  minimal twist d=%o;  good at 11: %o\n", dtwist, GoodReductionAt11(Cmin);
    catch e;
        printf "  MinimalTwist: %o\n", e`Object;
    end try;

    // Conductor exponents at odd primes up to TwistBound.
    badodd := [];
    for ell in PrimesUpTo(TwistBound) do
        if ell eq 2 then continue; end if;
        for p in PrimesAbove(K, ell) do
            e := ConductorExponentAt(Cmin, p);
            if e gt 0 then
                pnm := Type(K) eq FldRat select ell else Norm(p);
                Append(~badodd, <pnm, e>);
            end if;
        end for;
    end for;
    printf "  bad odd primes (norm:exp, up to %o): %o\n", TwistBound, badodd;

    // Frobenius polynomial at 11 and its mod-11 factorization.
    for P11 in PrimesAbove(K, 11) do
        try
            fp  := FrobPolyAt(Cmin, P11);
            fp11 := R11 ! fp;
            rts  := [r[1] : r in Roots(fp11)];
            pnm  := Type(K) eq FldRat select 11 else Norm(P11);
            printf "  Frob at P|11 (norm %o): %o = %o  (mod-11 roots in F_11: %o)\n",
                pnm, fp, Factorization(fp), rts;
        catch e;
            printf "  Frob at 11 not computable: %o\n", e`Object;
        end try;
    end for;
    printf "\n";
end for;
