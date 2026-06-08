// search_gp11.m -- Screening loop for the (1,11) Gross-Popescu pipeline.
//
// Enumerates rational base points a = (a1:a2:a3:a4) in P^3(Q) with |ai| <= BoxM,
// computes g(t) = Pf(S)(a1,a2,a3,a4,t) for each, factors g over Q, and for every
// irreducible factor h of degree <= DegBound that has a real root, tests the D_2
// point [a1:a2:a3:a4:t0] (t0 a root of h, over K = Q[t]/(h)) as an (1,11) surface.
//
// Filters (cheap, on the period matrix / Igusa invariants):
//   - non-degenerate D_2 point: Rank(S) = 4 (checked by PencilFromTorsion)
//   - inversion convergence: |R| below basin threshold
//   - potential good reduction at 11 (Igusa valuation criterion)
//   - geometric End = Z (simple abelian surface)
// Survivors stream to survivors_gp11.m for the heavy analysis (analyze_gp11.m).
//
// Run: magma search_gp11.m
AttachSpec("../CHIMP/CHIMP.spec");
import "GP11.m": IntertwiningMatrix, PencilFromTorsion;
import "InversionGP11.m": InvertGP11Fast;
import "Genus2Curve.m": IgusaInvariantsInK, GeometricEndomorphismDimension;
import "Reduction.m": PrimesAbove;

// ---- configuration ----
BoxM     := 2;    // |a_i| <= BoxM for i=1..4 (rational base point for P^3)
DegBound := 4;    // max degree of irreducible Pfaffian factors to consider
MaxPts   := 60;   // max P^3(Q) base points to screen
Prec     := 100;  // precision for inversion and invariant recognition
Trials   := 80;   // Levenberg-Marquardt restarts per D_2 point
survfile := "survivors_gp11.m";
logfile  := "search_results_gp11.txt";
// -----------------------

// Igusa valuation criterion for potential good reduction at 11.
function PotentialGoodReductionAt11(QI, K)
    wts := [2,4,6,8,10];
    for p in PrimesAbove(K, 11) do
        if QI[5] eq 0 then return false; end if;
        base := Valuation(QI[5], p) / 10;
        for i in [1..4] do
            if QI[i] ne 0 and Valuation(QI[i],p)/wts[i] lt base then return false; end if;
        end for;
    end for;
    return true;
end function;

// Norm-product of bad primes away from {2,11} up to Bound, as a proxy for conductor size.
function ConductorProxy11(QI, K : Bound := 100)
    wts := [2,4,6,8,10]; bad := []; sz := 1;
    for ell in PrimesUpTo(Bound) do
        if ell eq 2 or ell eq 11 then continue; end if;
        for p in PrimesAbove(K, ell) do
            if QI[5] eq 0 then continue; end if;
            base := Valuation(QI[5],p) / 10;
            if exists{i : i in [1..4] | QI[i] ne 0 and Valuation(QI[i],p)/wts[i] lt base} then
                Append(~bad, p);
                sz *:= Type(K) eq FldRat select ell else Norm(p);
            end if;
        end for;
    end for;
    return sz, bad;
end function;

// Primitive P^3(Q) points with max |coord| <= M, deduplicated and sorted by height.
function EnumerateP3Q(M)
    seen := {}; res := [];
    for a1 in [-M..M], a2 in [-M..M], a3 in [-M..M], a4 in [-M..M] do
        x := [Integers()| a1, a2, a3, a4];
        nz := [i : i in [1..4] | x[i] ne 0];
        if #nz eq 0 then continue; end if;
        g := GCD([Abs(c) : c in x]);
        x := [c div g : c in x];
        j := nz[#nz]; if x[j] lt 0 then x := [-c : c in x]; end if;
        key := <x[1], x[2], x[3], x[4]>;
        if key in seen then continue; end if;
        Include(~seen, key);
        Append(~res, <x, Max([Abs(c) : c in x])>);
    end for;
    Sort(~res, func<a,b | a[2]-b[2]>);
    return [r[1] : r in res];
end function;

Qt<t> := PolynomialRing(Rationals());
pts := EnumerateP3Q(BoxM);
printf "Enumerated %o P^3(Q) base points (BoxM=%o); testing up to %o.\n",
    #pts, BoxM, Min(MaxPts, #pts);
System("rm -f " cat survfile cat " " cat logfile);
PrintFile(survfile,
    "// Each entry: < [a1,a2,a3,a4], h_coeffs, QI_coordseqs, height >\n" cat
    "// h_coeffs: Q-coefficients of min-poly of x5, degree 0 upward.\n" cat
    "// QI_coordseqs: 5 lists, each the K-coordinate sequence of one Igusa invariant.\n" cat
    "survivors_gp11 := [* *];");
PrintFile(logfile, Sprintf("# (1,11) GP search: BoxM=%o DegBound=%o MaxPts=%o Prec=%o",
    BoxM, DegBound, MaxPts, Prec));
nsurv := 0;

for idx in [1..Min(MaxPts, #pts)] do
    a := pts[idx]; ht := Max([Abs(c) : c in a]);
    g := Pfaffian(IntertwiningMatrix([Qt| a[1],a[2],a[3],a[4],t]));
    if g eq 0 then continue; end if;
    facts := [f[1] : f in Factorization(g) | Degree(f[1]) ge 1 and Degree(f[1]) le DegBound];
    if #facts eq 0 then continue; end if;
    printf "[%o/%o] a=%o  (ht=%o, %o factor(s) of deg <= %o)\n",
        idx, Min(MaxPts,#pts), a, ht, #facts, DegBound;

    for h in facts do
        d := Degree(h);
        rts_real := [r[1] : r in Roots(h, RealField(Prec))];
        if #rts_real eq 0 then printf "  h=deg%o: no real roots, skip\n", d; continue; end if;

        CC := ComplexFieldExtra(Prec);
        if d eq 1 then
            K := Rationals();
            t0 := -Coefficient(h,0) / LeadingCoefficient(h);
            Ptors := [K| a[1], a[2], a[3], a[4], t0];
            emb := func<x | CC!x>;
        else
            K<aa> := NumberField(h);
            rt := CC ! rts_real[1];
            emb := hom<K -> CC | rt>;
            Ptors := [K| a[1], a[2], a[3], a[4], aa];
        end if;

        printf "  h=deg%o  K=%o\n", d, Type(K) eq FldRat select "Q" else Sprint(DefiningPolynomial(K));

        // degeneracy check: PencilFromTorsion errors if Rank(S) != 4
        v := []; w := []; ok := true;
        try v, w := PencilFromTorsion(Ptors);
        catch e; printf "    degenerate (%o)\n", e`Object; ok := false; end try;
        if not ok then continue; end if;

        // inversion
        tau := 0; QI := []; ok := true;
        try
            vCC := [emb(c) : c in v]; wCC := [emb(c) : c in w];
            tau, nr := InvertGP11Fast(vCC, wCC, CC : trials := Trials);
            printf "    |R| = %o\n", RealField(6)!nr;
            QI := IgusaInvariantsInK(tau, K, emb);
        catch e; printf "    inversion/recognition failed: %o\n", e`Object; ok := false; end try;
        if not ok then continue; end if;

        // filter 1: potential good reduction at 11
        if not PotentialGoodReductionAt11(QI, K) then
            printf "    potential-bad at 11, skip\n"; continue;
        end if;

        // filter 2: geometric End = Z
        ed := -1; try ed := GeometricEndomorphismDimension(tau); catch e; end try;
        if ed ne 1 then printf "    End-dim=%o (not simple), skip\n", ed; continue; end if;

        // field of moduli: check if invariants actually land in Q
        if Type(K) ne FldRat and forall{q : q in QI | Degree(MinimalPolynomial(q)) eq 1} then
            fom_str := "Q (descended)";
        else
            fom_str := Type(K) eq FldRat select "Q" else "K=" cat Sprint(DefiningPolynomial(K));
        end if;

        sz, bad := ConductorProxy11(QI, K);
        nsurv +:= 1;
        printf "    SURVIVOR  fom=%o  conductor-proxy=%o\n", fom_str, sz;

        // serialize invariants as coordinate-sequences over Q
        if d eq 1 then
            QI_seqs := [[Rationals()!q] : q in QI];
        else
            QI_seqs := [[Rationals()!c : c in Eltseq(q)] : q in QI];
        end if;
        h_coeffs := Coefficients(h);   // from degree 0 upward
        PrintFile(survfile, Sprintf(
            "Append(~survivors_gp11, < %o, %o, %o, %o >);",
            a, h_coeffs, QI_seqs, ht));
        PrintFile(logfile, Sprintf(
            "SURVIVOR a=%o h=%o ht=%o fom=%o proxy=%o bad=%o",
            a, h, ht, fom_str, sz, bad));
    end for;
end for;

printf "=== %o survivors recorded to %o; run analyze_gp11.m for the full analysis ===\n",
    nsurv, survfile;
PrintFile(logfile, Sprintf("# done: %o survivors", nsurv));
