// search_quadratic_ext.m -- Extends the Q(sqrt(2)) (1,5) HM search.
//
// Continues from where search_quadratic.m stopped:
//   - Re-runs points 1..30 with Trials=300 (to catch prior convergence failures).
//   - Then screens points 31..MaxPts with Trials=200.
//   - Deduplicates by Igusa invariants against existing survivors.m entries.
//   - Appends only genuinely new survivors to survivors.m.
//
// Run: magma search_quadratic_ext.m

AttachSpec("../CHIMP/CHIMP.spec");
load "survivors.m";   // defines k<s2> and survivors
import "InversionFast.m": x_to_tau_fast;
import "Genus2Curve.m": IgusaInvariantsInK, GeometricEndomorphismDimension;
import "Reduction.m": PotentialGoodReductionAt5, ConductorSizeProxy, PrimesAbove;
import "Heights.m": BoundedHeightPoints;

// ---- configuration ----
BoxM     := 2;
MaxPts   := 300;   // screen this many total non-rational K-points
RetryEnd := 30;    // re-screen these (with extra trials) to catch past failures
TrialsRetry := 300;
TrialsNew   := 200;
Prec     := 200;
survfile := "survivors.m";
logfile  := "search_results.txt";
// -----------------------

// Build set of already-found Igusa invariant tuples for deduplication.
seen_QI := [ s[2] : s in survivors ];

function IsKnown(QI, seen)
    return exists{ sq : sq in seen | forall{ j : j in [1..5] | QI[j] eq sq[j] } };
end function;

pts, hts := BoundedHeightPoints(k, BoxM : ExcludeRational := true);
total := Min(MaxPts, #pts);
printf "Enumerated %o non-rational K-points (BoxM=%o); screening up to %o.\n", #pts, BoxM, total;

nsurv := 0;
for i in [1..total] do
    x    := pts[i];
    tr   := i le RetryEnd select TrialsRetry else TrialsNew;
    printf "[%o/%o] h=%o x=%o trials=%o ... ", i, total, hts[i], x, tr;
    tau := 0; QI := []; ok := true;
    try
        CC := ComplexFieldExtra(Prec);
        tau, _, emb := x_to_tau_fast(x, CC : K := k, trials := tr);
        QI := IgusaInvariantsInK(tau, k, emb);
    catch e; printf "skip (%o)\n", e`Object; ok := false; end try;
    if not ok then continue; end if;
    if IsKnown(QI, seen_QI) then printf "already found\n"; continue; end if;
    if forall{ q : q in QI | Eltseq(k!q)[2] eq 0 } then printf "over Q (base change)\n"; continue; end if;
    if not PotentialGoodReductionAt5(QI, k) then printf "potential-bad at 5\n"; continue; end if;
    ed := -1; try ed := GeometricEndomorphismDimension(tau); catch e; end try;
    if ed ne 1 then printf "End-dim %o (not simple)\n", ed; continue; end if;
    size, bad := ConductorSizeProxy(QI, k);
    nsurv +:= 1;
    printf "NEW SURVIVOR (genuine, simple, pot-good@5); conductor-proxy=%o\n", size;
    PrintFile(survfile, Sprintf(
        "Append(~survivors, < [k| %o, %o, %o, %o], [k| %o, %o, %o, %o, %o], %o >);",
        x[1],x[2],x[3],x[4], QI[1],QI[2],QI[3],QI[4],QI[5], hts[i]));
    PrintFile(logfile, Sprintf(
        "SURVIVOR x=%o h=%o conductor-proxy=%o [EXT i=%o]", x, hts[i], size, i));
    Append(~seen_QI, QI);
end for;
printf "=== %o new survivors appended to %o ===\n", nsurv, survfile;
PrintFile(logfile, Sprintf("# ext done: %o new survivors", nsurv));
