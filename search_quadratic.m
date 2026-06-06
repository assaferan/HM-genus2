// Search driver: over a family of points x in P^3(K), screen for genus-2 curves with
// GOOD REDUCTION AT 5 and rank the survivors by a cheap conductor-size proxy.
//
// Every curve from this family already has a rational cyclic 5-isogeny, so the search
// is over good reduction at 5 (the hard filter) and small conductor.
//
// Per-x cost is the period computation for the Igusa invariants (Genus2Invariants:
// the stochastic basin search, no Mestre descent). The good-reduction filter and the
// conductor-size proxy are then cheap and model-free. Survivors are ranked; for the
// best few you can build the model with Genus2Curve and read off exact conductor
// exponents at odd primes with Reduction.m's ReductionReport.
//
// Run:  magma search_quadratic.m       (long: minutes per x)
// Customise K, the family `xs`, Prec and Trials below.

AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2Invariants;
import "Reduction.m": PotentialGoodReductionAt5, ConductorSizeProxy, PrimesAbove;

// ---------------- configuration ----------------
k<s2> := QuadraticField(2);
Prec   := 200;
Trials := 120;
// family of x in P^3(k) to screen (first coordinate fixed to 1 to pin the scaling):
xs := [ [k| 1, s2, c, d] : c, d in [2,3,4,5] ];
outfile := "search_results.txt";
// ------------------------------------------------

System("rm -f " cat outfile);
PrintFile(outfile, Sprintf("# good-reduction-at-5 search over %o, primes above 5 have norm %o",
          k, [Norm(p) : p in PrimesAbove(k,5)]));

survivors := [* *];
for i in [1..#xs] do
    x := xs[i];
    printf "[%o/%o] x = %o ... ", i, #xs, x;
    ok := true; QI := [];
    try
        QI := Genus2Invariants(x : K := k, Prec := Prec, Trials := Trials);
    catch e
        printf "invariants FAILED (%o)\n", e`Object; ok := false;
    end try;
    if not ok then continue; end if;
    if not PotentialGoodReductionAt5(QI, k) then
        printf "bad at 5 -> skip\n"; continue;
    end if;
    size, bad := ConductorSizeProxy(QI, k);
    printf "GOOD at 5; conductor-size(away 2,5)=%o badprime-norms=%o\n", size, [Norm(p) : p in bad];
    Append(~survivors, <x, QI, size, [Norm(p) : p in bad]>);
    PrintFile(outfile, Sprintf("GOOD5 x=%o size=%o badnorms=%o QI=%o", x, size, [Norm(p):p in bad], QI));
end for;

// rank survivors by conductor-size proxy (smaller is nicer)
surv := [ s : s in survivors ];
Sort(~surv, func< a, b | a[3] - b[3] >);
printf "\n=== %o survivor(s) with good reduction at 5, ranked by conductor-size (away from 2,5) ===\n", #surv;
for s in surv do
    printf "x = %o   size = %o   bad-prime norms = %o\n", s[1], s[3], s[4];
end for;
PrintFile(outfile, Sprintf("# done: %o survivors", #surv));
