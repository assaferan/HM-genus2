// Search driver (SCREENING): enumerate non-rational K-points x of P^3 in increasing
// height order and keep those whose genus-2 curve over K = Q(sqrt2) is
//   - genuinely over K (not a base change from Q: x not in P^3(Q) AND Igusa invariants
//     not all rational),
//   - of potential good reduction at 5, and
//   - SIMPLE (End(J) = Z).
// These three filters are cheap (period computation + invariant/endomorphism checks);
// the EXPENSIVE per-survivor analysis (Mestre descent, minimal twist, conductor,
// Frobenius-at-5) is deferred to analyze_survivors.m, which reuses the recognized
// invariants saved here (no second period computation).
//
// Survivors stream to survivors.m (loadable: k and a list <x, Igusa invariants, height>)
// and a human log to search_results.txt.  Run:  magma search_quadratic.m

AttachSpec("../CHIMP/CHIMP.spec");
import "InversionFast.m": x_to_tau_fast;
import "Genus2Curve.m": IgusaInvariantsInK, GeometricEndomorphismDimension;
import "Reduction.m": PotentialGoodReductionAt5, ConductorSizeProxy, PrimesAbove;
import "Heights.m": BoundedHeightPoints;

// ---------------- configuration ----------------
k<s2>   := QuadraticField(2);
BoxM    := 2;     // coordinate box for the height enumeration
MaxPts  := 30;    // how many lowest-height non-rational points to screen
Prec    := 200;   // precision (enough to recognize the invariants; descent is deferred)
Trials  := 120;   // basin-search restarts per point
survfile := "survivors.m";
logfile  := "search_results.txt";
// ------------------------------------------------

pts, hts := BoundedHeightPoints(k, BoxM : ExcludeRational := true);
printf "Enumerated %o non-rational points (BoxM=%o); screening %o lowest by height.\n", #pts, BoxM, Min(MaxPts,#pts);
System("rm -f " cat survfile cat " " cat logfile);
PrintFile(survfile, "k<s2> := QuadraticField(2); survivors := [* *];");
PrintFile(logfile, Sprintf("# screening over %o: genuine-K, potential-good-at-5, simple (End=Z)", k));

nsurv := 0;
for i in [1 .. Min(MaxPts, #pts)] do
    x := pts[i];
    printf "[%o/%o] h=%o x=%o ... ", i, Min(MaxPts,#pts), hts[i], x;
    tau := 0; QI := []; ok := true;
    try
        CC := ComplexFieldExtra(Prec);
        tau, _, emb := x_to_tau_fast(x, CC : K := k, trials := Trials);
        QI := IgusaInvariantsInK(tau, k, emb);
    catch e; printf "skip (%o)\n", e`Object; ok := false; end try;
    if not ok then continue; end if;
    if forall{ q : q in QI | Eltseq(k!q)[2] eq 0 } then printf "defined over Q (base change)\n"; continue; end if;
    if not PotentialGoodReductionAt5(QI, k) then printf "potential-bad at 5\n"; continue; end if;
    ed := -1; try ed := GeometricEndomorphismDimension(tau); catch e; end try;
    if ed ne 1 then printf "End-dim %o (not simple)\n", ed; continue; end if;
    size, bad := ConductorSizeProxy(QI, k);
    nsurv +:= 1;
    printf "SURVIVOR (genuine, simple, pot-good@5); conductor-size proxy=%o\n", size;
    PrintFile(survfile, Sprintf("Append(~survivors, < [k| %o, %o, %o, %o], [k| %o, %o, %o, %o, %o], %o >);",
        x[1],x[2],x[3],x[4], QI[1],QI[2],QI[3],QI[4],QI[5], hts[i]));
    PrintFile(logfile, Sprintf("SURVIVOR x=%o h=%o conductor-size-proxy=%o", x, hts[i], size));
end for;
printf "=== %o survivors recorded to %o (analyze with analyze_survivors.m) ===\n", nsurv, survfile;
PrintFile(logfile, Sprintf("# done: %o survivors", nsurv));
