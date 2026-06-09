// analyze_rational.m -- Per-survivor analysis for the (1,5) rational search.
// Reads survivors_rational.m (produced by search_rational.m), deduplicates by Igusa
// invariants, then for each distinct curve: builds the curve over Q, applies the
// minimal twist, reports the true conductor, and prints the Frobenius polynomial at 5.
//
// Run:  magma analyze_rational.m

AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, FrobeniusPolynomial, ConductorExponentAt,
                      GoodReductionAt5, PrimesAbove, ConductorSizeProxy;

load "survivors_rational.m";    // defines survivors_rational := [* <x, QI, ht>, ... *]
K := Rationals();
TwistBound := 50;
ZT<T> := PolynomialRing(Integers());
F5 := GF(5); R5<T5> := PolynomialRing(F5);

// Deduplicate by WPS-normalized Igusa invariants.
seen := {}; uniq := [* *];
for s in survivors_rational do
    key := < s[2][1], s[2][2], s[2][3], s[2][4], s[2][5] >;
    if key notin seen then
        Include(~seen, key);
        Append(~uniq, s);
    end if;
end for;

printf "survivors_rational: %o total, %o distinct by invariants.\n",
    #survivors_rational, #uniq;

for s in uniq do
    x := s[1]; QI := s[2]; ht := s[3];
    size, bad := ConductorSizeProxy(QI, K);
    printf "\n=== x = %o  (ht=%o, proxy=%o, bad=%o) ===\n", x, ht, size, bad;
    printf "  Igusa invariants: %o\n", QI;

    ok := true;
    try C := Genus2CurveFromIgusa(QI, K);
    catch e; printf "  build failed: %o\n", e`Object; ok := false; end try;
    if not ok then continue; end if;
    printf "  C = %o\n", C;

    Cmin := C; d := 1;
    try
        Cmin, d := MinimalTwist(C, K : Bound := TwistBound);
        printf "  minimal twist d = %o\n", d;
    catch e; printf "  MinimalTwist: %o\n", e`Object; end try;

    // Conductor at all odd primes up to TwistBound.
    badodd := [];
    for ell in PrimesUpTo(TwistBound) do
        if ell eq 2 then continue; end if;
        e := ConductorExponentAt(Cmin, ell);
        if e gt 0 then Append(~badodd, <ell, e>); end if;
    end for;
    printf "  good at 5: %o   bad odd primes (p, exp): %o\n", GoodReductionAt5(Cmin), badodd;

    // Full conductor (works over Q).
    try
        N := Conductor(Cmin);
        printf "  conductor N = %o = %o\n", N, Factorization(N);
    catch e; printf "  Conductor: %o\n", e`Object; end try;

    // Frobenius at 5.
    try
        fp := FrobeniusPolynomial(Cmin, 5);
        printf "  Frob at 5: %o = %o  (mod-5 roots: %o)\n",
            fp, Factorization(fp), [r[1] : r in Roots(R5!fp)];
    catch e; printf "  Frob at 5: %o\n", e`Object; end try;
end for;

printf "\nDone.\n";
