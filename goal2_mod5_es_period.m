// goal2_mod5_es_period.m -- produce B_f's period matrix via the ES/Oda method, in-process,
// with the fixes needed for our form:
//   (1) bypass GuessConductor by passing KnownConductor = disc(F)^2 * Norm(level) * Norm(cond chi)^2
//       (valid for chi coprime to the level; we use ONLY coprime-conductor characters);
//   (2) per Hodge-sign class, pick the first character with a NON-VANISHING central L-value.
// Then assemble Omegas and reconstruct the small period matrix tau in H_2.
//
// Args: B (char search bound, default 60), MAXN (eigenvalue norm bound, default 10103).
// Run: magma B:=60 MAXN:=10103 goal2_mod5_es_period.m
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
AttachSpec("../EichlerShimuraHMF/src/spec");
import "Genus2Curve.m": IgusaNumericFromTau;
if assigned B then B := StringToInteger(B); else B := 60; end if;
if assigned MAXN then MAXN := StringToInteger(MAXN); else MAXN := 10103; end if;

label := "2.2.8.1-5000.1-j";
DIR := "/scratch/home/assaferan/GitHub/HM-genus2/es_eig/";
LOGF := "goal2_mod5_es_period_out.txt";
PrintFile(LOGF, Sprintf("# ES period matrix, B=%o maxn=%o", B, MAXN) : Overwrite:=true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

f := LMFDBHMFwithEigenvalues(label, DIR);
Mf := Parent(f); H := HeckeEigenvalueField(Mf); dim := Dimension(Mf); FQ := BaseField(Mf);
levelnorm := Norm(Level(Mf));
discF := Abs(Discriminant(Integers(FQ)));
cond0 := discF^2 * levelnorm;      // untwisted L-conductor = 64*5000 = 320000
LOG(Sprintf("disc(F)=%o, level norm=%o, untwisted L-conductor=%o", discF, levelnorm, cond0));

chis := QuadraticCharactersUpTo(B, FQ);
signs := [HodgeSigns(chi) : chi in chis];
LOG(Sprintf("%o characters up to B=%o", #chis, B));

// L-value with KnownConductor (coprime-conductor chars only)
function Lval(chi, emb)
    c := Norm(Conductor(chi));
    if Gcd(c, 10) ne 1 then return false, 0, 0; end if;     // skip non-coprime (formula invalid)
    kc := cond0 * c^2;
    try
        L := LSeriesTwisted(f, chi : Embedding:=emb, maxn:=MAXN, KnownConductor:=kc);
        LSetPrecision(L, MaxPrecision(L, MAXN));
        return true, Evaluate(L, 1), CFENew(L);
    catch err
        return false, 0, 0;
    end try;
end function;

// for each sign class, find the first coprime character with non-vanishing L(f,chi)(1)
desired := [ [1,1], [1,-1], [-1,1], [-1,-1] ];
res2 := [* *];
ok := true;
for s in desired do
    cand := [j : j in [1..#chis] | signs[j] eq s];
    sign_res := [* *];
    for j in cand do
        okc, sp1, e1 := Lval(chis[j], 1);
        cnorm := Norm(Conductor(chis[j]));
        if not okc then LOG(Sprintf("    sign %o char %o (cond %o): skipped/errored", s, j, cnorm)); continue; end if;
        prec := Precision(sp1);
        LOG(Sprintf("    sign %o char %o (cond %o): |L1|=%o (prec %o)", s, j, cnorm, RealField(6)!Abs(sp1), prec));
        if Abs(sp1) le 10^-(prec*0.6) then continue; end if;     // vanishing -> skip
        // compute the other embedding too
        _, sp2, e2 := Lval(chis[j], 2);
        omset := [* *]; errset := [* *];
        for k->sp in [* sp1, sp2 *] do
            CC := ComplexFieldExtra(Precision(sp));
            om := -4*Pi(CC)^2*Sqrt(CC!2)*GaussSumOdaSimpleModuloSign(chis[j], s, CC)*(CC!sp);
            Append(~omset, om); Append(~errset, Abs(k eq 1 select e1 else e2));
        end for;
        Append(~sign_res, <j, <o : o in omset>, <er : er in errset> >);
        LOG(Sprintf("  sign %o: using char %o (cond %o), |L1|=%o", s, j, Norm(Conductor(chis[j])), RealField(6)!Abs(sp1)));
        break;
    end for;
    if #sign_res eq 0 then LOG(Sprintf("  sign %o: NO usable coprime non-vanishing character (raise B)", s)); ok := false; end if;
    Append(~res2, sign_res);
end for;
if not ok then LOG("cannot form all 4 Omega components; aborting"); exit; end if;

Omegas_per_sign := [* [* elt[2] : elt in vps *] : vps in res2 *];
Omegas := OmegasViaCremonaTrick(H, Omegas_per_sign);
raw := PossibleModuliPoints(H, Omegas);
possible_zs := [[ FixModuliPoint(H, z) : z in zs] : zs in raw];
LOG(Sprintf("possible moduli-point groups: %o", #possible_zs));

OH := Integers(H);
for gi->zs in possible_zs do
    for zi->z in zs do
        try tau := SmallPeriodMatrix(z, 1*OH, 1*OH);
        catch err LOG(Sprintf("  z[%o,%o]: SmallPeriodMatrix failed", gi, zi)); continue; end try;
        LOG(Sprintf("=== PERIOD MATRIX  z[%o,%o] ===", gi, zi));
        LOG(Sprintf("%o", ChangeRing(tau, ComplexField(10))));
        try IC := IgusaNumericFromTau(tau);
            LOG(Sprintf("  Igusa (numeric) = %o", [ComplexField(10)!x : x in IC]));
        catch err LOG(Sprintf("  IgusaNumericFromTau failed (precision): %o", err`Object)); end try;
    end for;
end for;
LOG("done");
exit;
