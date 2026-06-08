// search_rational.m -- (1,5) screening over Q.
//
// Analogue of search_quadratic.m for the rational case: enumerate P^3(Q) points x of
// bounded height, invert x -> tau -> Igusa invariants (over Q), and keep those with
// potential good reduction at 5 and geometric End = Z.
// No "genuinely over K" filter: all points are rational, so all curves are over Q.
//
// Run: magma search_rational.m
AttachSpec("../CHIMP/CHIMP.spec");
import "InversionFast.m": x_to_tau_fast;
import "Genus2Curve.m": IgusaInvariantsInK, GeometricEndomorphismDimension;
import "Reduction.m": PotentialGoodReductionAt5, ConductorSizeProxy;

// ---- configuration ----
BoxM     := 2;    // |a_i| <= BoxM for i=1..4
MaxPts   := 50;   // max points to screen
Prec     := 200;  // precision (200 is enough to recognize rational Igusa invariants)
Trials   := 80;   // random-restart trials per point
survfile := "survivors_rational.m";
logfile  := "search_results_rational.txt";
// -----------------------

// Primitive P^3(Q) points with max |coord| <= M, sorted by height.
function EnumerateP3Q(M)
    seen := {}; res := [];
    for a1 in [-M..M], a2 in [-M..M], a3 in [-M..M], a4 in [-M..M] do
        x := [Integers()| a1,a2,a3,a4];
        nz := [i : i in [1..4] | x[i] ne 0];
        if #nz eq 0 then continue; end if;
        g := GCD([Abs(c) : c in x]);
        x := [c div g : c in x];
        j := nz[#nz]; if x[j] lt 0 then x := [-c : c in x]; end if;
        key := <x[1],x[2],x[3],x[4]>;
        if key in seen then continue; end if;
        Include(~seen, key);
        Append(~res, <x, Max([Abs(c) : c in x])>);
    end for;
    Sort(~res, func<a,b | a[2]-b[2]>);
    return [r[1] : r in res], [r[2] : r in res];
end function;

K := Rationals();
pts, hts := EnumerateP3Q(BoxM);
printf "Enumerated %o P^3(Q) points (BoxM=%o); screening %o.\n",
    #pts, BoxM, Min(MaxPts, #pts);
System("rm -f " cat survfile cat " " cat logfile);
PrintFile(survfile, "survivors_rational := [* *];");
PrintFile(logfile, Sprintf("# (1,5) rational search: pot-good@5, simple (End=Z)"));
nsurv := 0;

for i in [1..Min(MaxPts, #pts)] do
    x := pts[i]; ht := hts[i];
    printf "[%o/%o] h=%o x=%o ... ", i, Min(MaxPts,#pts), ht, x;
    tau := 0; QI := []; emb := 0; ok := true;
    try
        CC := ComplexFieldExtra(Prec);
        tau, _, emb := x_to_tau_fast([K|c : c in x], CC : K := K, trials := Trials);
        QI := IgusaInvariantsInK(tau, K, emb);
    catch e; printf "skip (%o)\n", e`Object; ok := false; end try;
    if not ok then continue; end if;
    if not PotentialGoodReductionAt5(QI, K) then printf "potential-bad at 5\n"; continue; end if;
    ed := -1; try ed := GeometricEndomorphismDimension(tau); catch e; end try;
    if ed ne 1 then printf "End-dim %o (not simple)\n", ed; continue; end if;
    size, bad := ConductorSizeProxy(QI, K);
    nsurv +:= 1;
    printf "SURVIVOR (pot-good@5, simple); conductor-proxy=%o bad=%o\n", size, bad;
    PrintFile(survfile, Sprintf(
        "Append(~survivors_rational, < [%o,%o,%o,%o], [%o,%o,%o,%o,%o], %o >);",
        x[1],x[2],x[3],x[4], QI[1],QI[2],QI[3],QI[4],QI[5], ht));
    PrintFile(logfile, Sprintf("SURVIVOR x=%o h=%o proxy=%o bad=%o", x, ht, size, bad));
end for;

printf "=== %o survivors recorded to %o ===\n", nsurv, survfile;
PrintFile(logfile, Sprintf("# done: %o survivors", nsurv));
