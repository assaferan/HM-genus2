// process_rational_d2.m -- For each rational non-degenerate D_2 point, run the (1,11)
// pipeline to get a genus-2 curve over Q, then filter: field of moduli Q, good reduction
// at 11, End=Z, and compute the conductor (over Q, Magma's Conductor works directly).
//
// Run: magma process_rational_d2.m   (needs rational_d2_points.txt from prescan)

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "GP11.m": PencilFromTorsion;
import "InversionGP11.m": InvertGP11Fast;
import "Genus2Curve.m": IgusaInvariantsInK, GeometricEndomorphismDimension, Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, PrimesAbove;

load "rational_d2_points.txt";   // d2pts := [ [a1..a4,x5], ... ]
Prec := 60; Trials := 24;   // low cap: convergent points break early; bad points abandon fast
d2pts := [ d2pts[i] : i in [1..Min(9, #d2pts)] ];
CC := ComplexFieldExtra(Prec);
emb := func<x | CC!x>;

plog := "ratproc_progress.txt";   // PrintFile flushes immediately -> watchable progress
PrintFile(plog, "# rational D_2 pipeline progress" : Overwrite := true);
procedure LOG(s) PrintFile(plog, s); end procedure;

printf "Processing %o rational D_2 points...\n\n", #d2pts;
survivors := [];
for i in [1..#d2pts] do
    P := d2pts[i];
    printf "=== [%o/%o] D_2 point %o ===\n", i, #d2pts, P;
    LOG(Sprintf("[%o/%o] D2=%o : start", i, #d2pts, P));
    t0clock := Cputime();
    Ptors := [Rationals()| c : c in P];

    ok := true; v := []; w := [];
    try v, w := PencilFromTorsion(Ptors);
    catch e; printf "  degenerate: %o\n", e`Object; ok := false; end try;
    if not ok then continue; end if;

    tau := 0; nr := 0;
    try
        tau, nr := InvertGP11Fast([emb(c):c in v], [emb(c):c in w], CC : trials := Trials);
        printf "  |R| = %o\n", RealField(6)!nr;
    catch e; printf "  inversion failed: %o\n", e`Object; ok := false; end try;
    LOG(Sprintf("    inversion |R|=%o  (%os)", ok select RealField(6)!nr else "FAILED", Round(Cputime(t0clock))));
    if not ok or nr gt 1e-20 then printf "  did not converge, skip\n"; LOG("    -> no converge, skip"); continue; end if;

    QI := [];
    try QI := IgusaInvariantsInK(tau, Rationals(), emb);
    catch e; printf "  invariant recognition failed: %o\n", e`Object; continue; end try;
    // field of moduli Q?
    if not forall{q : q in QI | q in Rationals()} then
        printf "  invariants not rational -> field of moduli != Q, skip\n";
        LOG("    -> invariants not rational (fom != Q), skip"); continue;
    end if;
    printf "  field of moduli Q. Igusa (normalized): %o\n", QI;

    // good reduction at 11 (valuation criterion)
    goodat11 := true;
    if QI[5] eq 0 then goodat11 := false; else
        base := Valuation(Rationals()!QI[5], 11)/10;
        for j in [1..4] do
            if QI[j] ne 0 and Valuation(Rationals()!QI[j], 11)/[2,4,6,8][j] lt base then goodat11 := false; end if;
        end for;
    end if;
    printf "  potential good reduction at 11: %o\n", goodat11;

    // End = Z ?
    ed := -1; try ed := GeometricEndomorphismDimension(tau); catch e; end try;
    printf "  geometric End dimension: %o %o\n", ed, ed eq 1 select "(End=Z, simple)" else "(NOT simple)";

    // reconstruct curve over Q + conductor
    cond := 0; condstr := "?";
    try
        C := Genus2CurveFromIgusa(QI, Rationals());
        Cmin := MinimalTwist(C, Rationals());
        cond := Conductor(Cmin);
        condstr := Sprintf("%o = %o", cond, Factorization(cond));
    catch e; condstr := "conductor failed: " cat Sprint(e`Object); end try;
    printf "  conductor: %o\n", condstr;

    Append(~survivors, <P, QI, goodat11, ed, condstr>);
    LOG(Sprintf("    SURVIVOR fom=Q good@11=%o End-dim=%o cond=%o", goodat11, ed, condstr));
    printf "\n";
end for;

printf "=== Summary: %o curves over Q ===\n", #survivors;
for s in survivors do
    printf "  D2=%o  good@11=%o  End-dim=%o  cond=%o\n", s[1], s[3], s[4], s[5];
end for;
