// Heavy per-survivor analysis, decoupled from the (fast) screening in search_quadratic.m.
// Reads survivors.m (the recognized Igusa invariants saved by the search -- so NO period
// computation is repeated), and for each survivor: builds the curve over K from its
// invariants, takes the MINIMAL quadratic twist (good reduction at every prime of
// potential good reduction), and reports the conductor exponents at odd primes and the
// factored Frobenius polynomial at 5.
//
// Run:  magma analyze_survivors.m

AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, FrobeniusPolynomial, ConductorExponentAt, GoodReductionAt5, PrimesAbove;

load "survivors.m";          // defines k<s2> and  survivors := [* <x, QI, height>, ... *]
TwistBound := 30;            // examine odd primes up to this for twist / conductor
ZT<T> := PolynomialRing(Integers());
p5 := PrimesAbove(k, 5)[1];

printf "Analyzing %o survivor(s).\n", #survivors;
for s in survivors do
    x := s[1]; QI := s[2]; h := s[3];
    printf "\n=== x = %o   (height %o) ===\n", x, h;
    ok := true;
    try
        C := Genus2CurveFromIgusa(QI, k);
    catch e; printf "  build failed (%o)\n", e`Object; ok := false; end try;
    if not ok then continue; end if;
    if not (BaseRing(C) cmpeq k) then printf "  does not descend to K\n"; continue; end if;
    Cmin, d := MinimalTwist(C, k : Bound := TwistBound);
    printf "  minimal twist d = %o ; good reduction at 5: %o\n", d, GoodReductionAt5(Cmin);
    badodd := [];
    for ell in PrimesUpTo(TwistBound) do
        if ell eq 2 then continue; end if;
        for p in PrimesAbove(k, ell) do
            e := ConductorExponentAt(Cmin, p);
            if e gt 0 then Append(~badodd, <Norm(p), e>); end if;
        end for;
    end for;
    printf "  bad odd primes (norm^exp, up to %o): %o\n", TwistBound, badodd;
    try printf "  Frobenius polynomial at 5: %o\n", Factorization(FrobeniusPolynomial(Cmin, p5));
    catch e; printf "  Frobenius-at-5 not computable (%o)\n", e`Object; end try;
end for;
