// search_qsqrt5_gp11_ext.m -- Extends the Q(sqrt(5)) (1,11) GP11 search.
//
// Continues from where search_qsqrt5_gp11.m stopped:
//   - Re-runs points 1..30 with Trials=200 (to catch prior convergence failures).
//   - Then screens points 31..MaxPts with Trials=150.
//   - Deduplicates by Igusa invariants against existing survivors_qsqrt5_gp11.m.
//   - Appends only genuinely new survivors.
//
// Run: magma search_qsqrt5_gp11_ext.m

AttachSpec("../CHIMP/CHIMP.spec");
load "survivors_qsqrt5_gp11.m";   // defines k<s5> and survivors_qsqrt5_gp11
import "GP11.m": IntertwiningMatrix, PencilFromTorsion;
import "InversionGP11.m": InvertGP11Fast;
import "Genus2Curve.m": IgusaInvariantsInK, GeometricEndomorphismDimension;
import "Reduction.m": PrimesAbove, ConductorProxy11;
import "Heights.m": BoundedHeightPoints;

// ---- configuration ----
BoxM        := 2;
MaxPts      := 200;
RetryEnd    := 30;
TrialsRetry := 200;
TrialsNew   := 150;
Prec        := 100;
survfile    := "survivors_qsqrt5_gp11.m";
logfile     := "search_results_qsqrt5_gp11.txt";
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


function RealEmbeddings(K, CC)
    rts := [r[1] : r in Roots(DefiningPolynomial(K), RealField(Precision(CC)))];
    return [hom<K -> CC | CC!rt> : rt in rts];
end function;

// Deduplication
seen_QI := [ s[2] : s in survivors_qsqrt5_gp11 ];
function IsKnown(QI, seen)
    return exists{ sq : sq in seen | forall{ j : j in [1..5] | QI[j] eq sq[j] } };
end function;

Kt<t> := PolynomialRing(k);
pts, hts := BoundedHeightPoints(k, BoxM : ExcludeRational := true);
total := Min(MaxPts, #pts);
printf "Enumerated %o non-rational K-points (BoxM=%o); screening up to %o.\n", #pts, BoxM, total;

nsurv := 0;
for i in [1..total] do
    x  := pts[i];
    tr := i le RetryEnd select TrialsRetry else TrialsNew;
    printf "[%o/%o] h=%o x=%o trials=%o ... ", i, total, hts[i], x, tr;

    g := Pfaffian(IntertwiningMatrix([Kt| x[1],x[2],x[3],x[4],t]));
    if g eq 0 then printf "Pfaffian=0\n"; continue; end if;
    rts_k := [r[1] : r in Roots(g, k)];
    if #rts_k eq 0 then printf "no K-root\n"; continue; end if;
    printf "%o K-root(s)\n", #rts_k;

    for t0 in rts_k do
        Ptors := [k| x[1], x[2], x[3], x[4], t0];
        v := []; w_vec := []; ok := true;
        try v, w_vec := PencilFromTorsion(Ptors);
        catch e; printf "  t0~%o: degenerate\n", t0; ok := false; end try;
        if not ok then continue; end if;

        CC := ComplexFieldExtra(Prec);
        embs := RealEmbeddings(k, CC);
        tau := 0; QI := []; found_emb := false;
        for emb in embs do
            try
                vCC := [emb(c) : c in v]; wCC := [emb(c) : c in w_vec];
                tau, nr := InvertGP11Fast(vCC, wCC, CC : trials := tr);
                QI := IgusaInvariantsInK(tau, k, emb);
                found_emb := true;
                printf "  t0~%o: |R|=%o\n", RealField(6)!Real(emb(t0)), RealField(6)!nr;
                break;
            catch e; end try;
        end for;
        if not found_emb then printf "  all embeddings failed\n"; continue; end if;

        if IsKnown(QI, seen_QI) then printf "  already found\n"; continue; end if;
        if forall{ q : q in QI | Eltseq(k!q)[2] eq 0 } then
            printf "  over Q (base change)\n"; continue; end if;
        if not PotentialGoodReductionAt11(QI, k) then
            printf "  potential-bad at 11\n"; continue; end if;
        ed := -1; try ed := GeometricEndomorphismDimension(tau); catch e; end try;
        if ed ne 1 then printf "  End-dim=%o (not simple)\n", ed; continue; end if;

        sz, bad := ConductorProxy11(QI, k);
        nsurv +:= 1;
        printf "  NEW SURVIVOR; conductor-proxy=%o\n", sz;
        PrintFile(survfile, Sprintf(
            "Append(~survivors_qsqrt5_gp11, < [k| %o, %o, %o, %o, %o], [k| %o, %o, %o, %o, %o], %o >);",
            x[1],x[2],x[3],x[4],t0, QI[1],QI[2],QI[3],QI[4],QI[5], hts[i]));
        PrintFile(logfile, Sprintf(
            "SURVIVOR x=%o t0=%o h=%o proxy=%o [EXT i=%o]", x, t0, hts[i], sz, i));
        Append(~seen_QI, QI);
    end for;
end for;
printf "=== %o new survivors appended to %o ===\n", nsurv, survfile;
PrintFile(logfile, Sprintf("# ext done: %o new survivors", nsurv));
