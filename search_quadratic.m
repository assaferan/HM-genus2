// Search driver: enumerate K-rational points x of P^3 in INCREASING HEIGHT order and
// screen for genus-2 curves with GOOD REDUCTION AT 5, ranking survivors by a conductor
// -size proxy. (Low height tends to give low conductor, so height order surfaces the
// best candidates first.)
//
// Every curve from this family already has a rational cyclic 5-isogeny, so the search
// is over good reduction at 5 + small conductor. Per-x cost is the period computation
// for the Igusa invariants (Genus2Invariants: stochastic search, no Mestre descent);
// the good-reduction filter and the conductor-size proxy are cheap and model-free.
//
// NOTE on scale: P^3(K) has ~B^8 points of height <= B, so an exhaustive height bound
// is infeasible beyond tiny B. We therefore enumerate within a coordinate box (|a|,|b|
// <= BoxM for coordinates a+b*sqrt(d)), sort by height, and process the lowest MaxPts
// of them. Increase BoxM/MaxPts for a wider search (cost grows fast).
//
// Run:  magma search_quadratic.m        (long: minutes per processed point)

AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2Invariants;
import "Reduction.m": PotentialGoodReductionAt5, ConductorSizeProxy, PrimesAbove;
import "Heights.m": BoundedHeightPoints;

// ---------------- configuration ----------------
k<s2>  := QuadraticField(2);
BoxM   := 2;        // coordinate box for the enumeration (a+b*s2, |a|,|b| <= BoxM)
MaxPts := 40;       // how many lowest-height points to actually process
Prec   := 200;
Trials := 120;
outfile := "search_results.txt";
// ------------------------------------------------

pts, hts := BoundedHeightPoints(k, BoxM);
printf "Enumerated %o points in the box (BoxM=%o); processing the %o lowest by height.\n",
       #pts, BoxM, Min(MaxPts, #pts);
System("rm -f " cat outfile);
PrintFile(outfile, Sprintf("# good-reduction-at-5 search over %o; primes above 5 have norm %o",
          k, [Norm(p) : p in PrimesAbove(k,5)]));

survivors := [];
for i in [1 .. Min(MaxPts, #pts)] do
    x := pts[i];
    printf "[%o/%o] height %o  x=%o ... ", i, Min(MaxPts,#pts), hts[i], x;
    QI := []; ok := true;
    try
        QI := Genus2Invariants(x : K := k, Prec := Prec, Trials := Trials);
    catch e
        printf "skip (%o)\n", e`Object; ok := false;
    end try;
    if not ok then continue; end if;
    if not PotentialGoodReductionAt5(QI, k) then printf "bad at 5\n"; continue; end if;
    size, bad := ConductorSizeProxy(QI, k);
    printf "GOOD at 5; cond-size(away 2,5)=%o badnorms=%o\n", size, [Norm(p) : p in bad];
    Append(~survivors, <x, hts[i], size, [Norm(p) : p in bad], QI>);
    PrintFile(outfile, Sprintf("GOOD5 height=%o x=%o size=%o badnorms=%o QI=%o",
              hts[i], x, size, [Norm(p):p in bad], QI));
end for;

Sort(~survivors, func< a, b | a[3] - b[3] >);   // by conductor-size proxy
printf "\n=== %o survivor(s) with good reduction at 5, ranked by conductor-size (away 2,5) ===\n", #survivors;
for s in survivors do
    printf "x=%o  height=%o  cond-size=%o  bad-prime-norms=%o\n", s[1], s[2], s[3], s[4];
end for;
PrintFile(outfile, Sprintf("# done: %o survivors", #survivors));
