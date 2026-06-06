// Search driver: enumerate K-rational points x of P^3 in INCREASING HEIGHT order and
// find genus-2 curves over K = Q(sqrt2) that are
//   (i)  SIMPLE (End(J) = Z, generic -- twist-invariant, from the period matrix),
//   (ii) of potential good reduction at 5, and after the MINIMAL TWIST genuinely good
//        at 5 (and at every other prime of potential good reduction), so the conductor
//        is as small as possible (supported on the genuinely-bad primes, e.g. 2),
// reporting each survivor's conductor exponents and factored Frobenius polynomial at 5.
// (Every curve in this family already has a rational cyclic 5-isogeny.)
//
// Cost: the period computation per x (stochastic) dominates; the End=Z and potential-
// good-at-5 filters are cheaper and are applied before the (expensive) Mestre descent,
// which runs only for survivors. Tune the config block; results stream to search_results.txt.
//
// Run:  magma search_quadratic.m

AttachSpec("../CHIMP/CHIMP.spec");
import "InversionFast.m": x_to_tau_fast;
import "Genus2Curve.m": IgusaInvariantsInK, Genus2CurveFromTau, GeometricEndomorphismDimension;
import "Reduction.m": PotentialGoodReductionAt5, MinimalTwist, FrobeniusPolynomial,
                      ConductorExponentAt, GoodReductionAt5, PrimesAbove;
import "Heights.m": BoundedHeightPoints;

// ---------------- configuration ----------------
k<s2>      := QuadraticField(2);
BoxM       := 2;     // coordinate box for the height enumeration
MaxPts     := 40;    // how many lowest-height points to process
Prec       := 300;   // working precision
Trials     := 120;   // basin-search restarts per point
TwistBound := 50;    // examine odd primes up to this for the minimal twist / conductor
outfile    := "search_results.txt";
// ------------------------------------------------

pts, hts := BoundedHeightPoints(k, BoxM : ExcludeRational := true);   // drop x in P^3(Q)
printf "Enumerated %o non-rational points (BoxM=%o); processing %o lowest by height.\n", #pts, BoxM, Min(MaxPts,#pts);
System("rm -f " cat outfile);
PrintFile(outfile, Sprintf("# simple + good-at-5 search over %o; p|5 norm %o", k, [Norm(p):p in PrimesAbove(k,5)]));
ZT<T> := PolynomialRing(Integers());
p5 := PrimesAbove(k, 5)[1];
survivors := [];

for i in [1 .. Min(MaxPts, #pts)] do
    x := pts[i];
    printf "[%o/%o] h=%o x=%o ... ", i, Min(MaxPts,#pts), hts[i], x;
    tau := 0; QI := []; emb := 0; ok := true;
    try
        CC := ComplexFieldExtra(Prec);
        tau, _, emb := x_to_tau_fast(x, CC : K := k, trials := Trials);
        QI := IgusaInvariantsInK(tau, k, emb);
    catch e; printf "skip (%o)\n", e`Object; ok := false; end try;
    if not ok then continue; end if;
    // reject curves defined over Q (base changes / their twists): rational Igusa invariants
    if forall{ q : q in QI | Eltseq(k!q)[2] eq 0 } then printf "defined over Q (base change)\n"; continue; end if;
    if not PotentialGoodReductionAt5(QI, k) then printf "potential-bad at 5\n"; continue; end if;
    endim := -1; try endim := GeometricEndomorphismDimension(tau); catch e; end try;
    if endim ne 1 then printf "End-dim %o (not simple)\n", endim; continue; end if;

    // survivor: build the curve, take the minimal twist, read off conductor + Frob_5
    okc := true;
    try
        C := Genus2CurveFromTau(tau, k, emb);
        if not (BaseRing(C) cmpeq k) then printf "does not descend to k\n"; continue; end if;
        Cmin, d := MinimalTwist(C, k : Bound := TwistBound);
    catch e; printf "curve/twist failed (%o)\n", e`Object; okc := false; end try;
    if not okc then continue; end if;

    badodd := [];
    for ell in PrimesUpTo(TwistBound) do
        if ell eq 2 then continue; end if;
        for p in PrimesAbove(k, ell) do
            ce := ConductorExponentAt(Cmin, p);
            if ce gt 0 then Append(~badodd, <Norm(p), ce>); end if;
        end for;
    end for;
    frob5 := "?";
    try frob5 := Sprint(Factorization(FrobeniusPolynomial(Cmin, p5))); catch e; end try;
    g5 := GoodReductionAt5(Cmin);
    oddnorm := &*([1] cat [ t[1]^t[2] : t in badodd ]);
    printf "SIMPLE; good@5=%o; bad odd primes (norm^exp)=%o; Frob_5=%o\n", g5, badodd, frob5;
    Append(~survivors, <x, hts[i], oddnorm, badodd, frob5>);
    PrintFile(outfile, Sprintf("x=%o h=%o good5=%o oddcond=%o badodd=%o Frob5=%o",
              x, hts[i], g5, oddnorm, badodd, frob5));
end for;

Sort(~survivors, func< a, b | a[3] - b[3] >);   // by odd-part conductor norm
printf "\n=== %o simple survivors, ranked by odd-part conductor norm ===\n", #survivors;
for s in survivors do
    printf "x=%o  h=%o  odd-cond-norm=%o  bad-odd=%o  Frob_5=%o\n", s[1], s[2], s[3], s[4], s[5];
end for;
PrintFile(outfile, Sprintf("# done: %o survivors", #survivors));
