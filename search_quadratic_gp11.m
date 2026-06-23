// search_quadratic_gp11.m -- Screening loop for the (1,11) pipeline over K = Q(sqrt(2)).
//
// Direct analogue of search_quadratic.m for (1,5). Enumerates non-rational K-points
// (x1:x2:x3:x4) in P^3(K) with bounded height, computes the Pfaffian polynomial
// g(t) = Pf(S)(x1,x2,x3,x4,t) in K[t], finds K-rational roots t0, and for each
// D_2 point [x1:x2:x3:x4:t0] in P^4(K) inverts to a genus-2 curve.
//
// Filters (cheap, from the period matrix / Igusa invariants):
//   - non-degenerate D_2: Rank(S) = 4 (checked by PencilFromTorsion)
//   - inversion convergence
//   - genuinely over K: Igusa invariants not all rational
//   - potential good reduction at 11 (Igusa valuation criterion)
//   - geometric End = Z (simple abelian surface)
// Survivors (with recognized Igusa invariants in K) stream to survivors_quadratic_gp11.m
// for the heavy analysis; run analyze_quadratic_gp11.m.
//
// Run: magma search_quadratic_gp11.m
AttachSpec("../CHIMP/CHIMP.spec");
import "GP11.m": IntertwiningMatrix, PencilFromTorsion;
import "InversionGP11.m": InvertGP11Fast;
import "Genus2Curve.m": IgusaInvariantsInK, GeometricEndomorphismDimension;
import "Reduction.m": PrimesAbove, ConductorProxy11;
import "Heights.m": BoundedHeightPoints;

// ---- configuration ----
k<s2>    := QuadraticField(2);
BoxM     := 2;    // coordinate box for height enumeration
MaxPts   := 30;   // max non-rational P^3(K) base points to screen
Prec     := 100;  // precision for inversion and recognition
Trials   := 80;   // LM restarts per (D_2 point, embedding) pair
survfile := "survivors_quadratic_gp11.m";
logfile  := "search_results_quadratic_gp11.txt";
// -----------------------

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

// Real embeddings of K into CC.
function RealEmbeddings(K, CC)
    rts := [r[1] : r in Roots(DefiningPolynomial(K), RealField(Precision(CC)))];
    return [hom<K -> CC | CC!rt> : rt in rts];
end function;

Kt<t> := PolynomialRing(k);
pts, hts := BoundedHeightPoints(k, BoxM : ExcludeRational := true);
printf "Enumerated %o non-rational K-points (BoxM=%o); screening %o.\n",
    #pts, BoxM, Min(MaxPts, #pts);
System("rm -f " cat survfile cat " " cat logfile);
PrintFile(survfile, "k<s2> := QuadraticField(2); survivors_quadratic_gp11 := [* *];");
PrintFile(logfile, Sprintf("# (1,11) GP search over %o: genuine-K, pot-good@11, simple", k));
nsurv := 0;

for i in [1..Min(MaxPts, #pts)] do
    x := pts[i];
    printf "[%o/%o] h=%o x=%o ... ", i, Min(MaxPts,#pts), hts[i], x;

    // Pfaffian in t over K; find K-rational roots.
    g := Pfaffian(IntertwiningMatrix([Kt| x[1],x[2],x[3],x[4],t]));
    if g eq 0 then printf "Pfaffian=0, skip\n"; continue; end if;
    rts_k := [r[1] : r in Roots(g, k)];
    if #rts_k eq 0 then printf "no K-rational root\n"; continue; end if;
    printf "%o K-root(s)\n", #rts_k;

    for t0 in rts_k do
        Ptors := [k| x[1], x[2], x[3], x[4], t0];

        // Degeneracy check.
        v := []; w := []; ok := true;
        try v, w := PencilFromTorsion(Ptors);
        catch e; printf "  t0 ~ %o: degenerate\n", t0; ok := false; end try;
        if not ok then continue; end if;

        // Try each real embedding of K; keep first converging one.
        CC := ComplexFieldExtra(Prec);
        embs := RealEmbeddings(k, CC);
        tau := 0; QI := []; found_emb := false;
        for emb in embs do
            try
                vCC := [emb(c) : c in v]; wCC := [emb(c) : c in w];
                tau, nr := InvertGP11Fast(vCC, wCC, CC : trials := Trials);
                QI := IgusaInvariantsInK(tau, k, emb);
                found_emb := true;
                printf "  t0 ~ %o: |R|=%o\n",
                    RealField(6)!Real(emb(t0)), RealField(6)!nr;
                break;
            catch e; end try;
        end for;
        if not found_emb then printf "  all embeddings failed\n"; continue; end if;

        // Filter 1: genuinely over K (not a base change from Q).
        if forall{ q : q in QI | Eltseq(k!q)[2] eq 0 } then
            printf "  defined over Q (base change)\n"; continue;
        end if;

        // Filter 2: potential good reduction at 11.
        if not PotentialGoodReductionAt11(QI, k) then
            printf "  potential-bad at 11\n"; continue;
        end if;

        // Filter 3: geometric End = Z.
        ed := -1; try ed := GeometricEndomorphismDimension(tau); catch e; end try;
        if ed ne 1 then printf "  End-dim=%o (not simple)\n", ed; continue; end if;

        sz, bad := ConductorProxy11(QI, k);
        nsurv +:= 1;
        printf "  SURVIVOR (genuine, simple, pot-good@11); conductor-proxy=%o\n", sz;
        PrintFile(survfile, Sprintf(
            "Append(~survivors_quadratic_gp11, < [k| %o, %o, %o, %o, %o], [k| %o, %o, %o, %o, %o], %o >);",
            x[1],x[2],x[3],x[4],t0, QI[1],QI[2],QI[3],QI[4],QI[5], hts[i]));
        PrintFile(logfile, Sprintf(
            "SURVIVOR x=%o t0=%o h=%o proxy=%o bad=%o", x, t0, hts[i], sz, bad));
    end for;
end for;

printf "=== %o survivors recorded to %o; run analyze_quadratic_gp11.m ===\n", nsurv, survfile;
PrintFile(logfile, Sprintf("# done: %o survivors", nsurv));
