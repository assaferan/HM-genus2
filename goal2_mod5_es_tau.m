// goal2_mod5_es_tau.m -- produce B_f's (small) period matrix from the THREE computable Oda
// period components Omega_{++}, Omega_{+-}, Omega_{-+}.  The moduli point that SmallPeriodMatrix
// consumes is z0 = [Omega_{-+}/Omega_{++} per Hecke embedding] (Omega_{--} only fixes a unit
// scaling via the cross-product, which we skip here => the matrix is correct up to that unit).
//
// Args: B (char bound, default 60), MAXN (default 10103).
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
AttachSpec("../EichlerShimuraHMF/src/spec");
import "Genus2Curve.m": IgusaNumericFromTau;
if assigned B then B := StringToInteger(B); else B := 60; end if;
if assigned MAXN then MAXN := StringToInteger(MAXN); else MAXN := 10103; end if;

label := "2.2.8.1-5000.1-j";
DIR := "/scratch/home/assaferan/GitHub/HM-genus2/es_eig/";
LOGF := "goal2_mod5_es_tau_out.txt";
PrintFile(LOGF, Sprintf("# raw ES period matrix from 3 components, B=%o maxn=%o", B, MAXN) : Overwrite:=true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

f := LMFDBHMFwithEigenvalues(label, DIR);
Mf := Parent(f); H := HeckeEigenvalueField(Mf); dim := Dimension(Mf); FQ := BaseField(Mf);
OH := Integers(H);
cond0 := Abs(Discriminant(Integers(FQ)))^2 * Norm(Level(Mf));   // 320000
chis := QuadraticCharactersUpTo(B, FQ); signs := [HodgeSigns(chi) : chi in chis];

function Lval(chi, emb)
    c := Norm(Conductor(chi));
    if Gcd(c, 10) ne 1 then return false, 0; end if;
    try
        L := LSeriesTwisted(f, chi : Embedding:=emb, maxn:=MAXN, KnownConductor:=cond0*c^2);
        LSetPrecision(L, MaxPrecision(L, MAXN));
        return true, Evaluate(L, 1);
    catch err return false, 0; end try;
end function;

// Omega for a sign class: first coprime non-vanishing char, both Hecke embeddings
function getOmega(s)
    for j in [k : k in [1..#chis] | signs[k] eq s] do
        okc, sp1 := Lval(chis[j], 1);
        if not okc then continue; end if;
        if Abs(sp1) le 10^-(Precision(sp1)*0.6) then continue; end if;
        _, sp2 := Lval(chis[j], 2);
        om := [];
        for sp in [sp1, sp2] do
            CC := ComplexFieldExtra(Precision(sp));
            Append(~om, -4*Pi(CC)^2*Sqrt(CC!2)*GaussSumOdaSimpleModuloSign(chis[j], s, CC)*(CC!sp));
        end for;
        return true, om, j;
    end for;
    return false, [], 0;
end function;

ok1, om_pp, j1 := getOmega([1,1]);
ok2, om_pm, j2 := getOmega([1,-1]);
ok3, om_mp, j3 := getOmega([-1,1]);
LOG(Sprintf("Omega components found: ++ (char %o), +- (char %o), -+ (char %o)", j1, j2, j3));
if not (ok1 and ok2 and ok3) then LOG("missing a component; abort"); exit; end if;

// raw moduli points (no Omega_-- unit rescaling)
prec := Min([Precision(x) : x in om_pp cat om_pm cat om_mp]);
z0 := [ om_mp[i]/om_pp[i] : i in [1..#om_pp] ];   // Omega_-+/Omega_++
z1 := [ om_pm[i]/om_pp[i] : i in [1..#om_pp] ];   // Omega_+-/Omega_++
LOG(Sprintf("working precision ~%o digits", prec));
LOG(Sprintf("z0 = Omega_-+/Omega_++ per embedding = %o", [ComplexField(6)!x : x in z0]));
LOG(Sprintf("z1 = Omega_+-/Omega_++ per embedding = %o", [ComplexField(6)!x : x in z1]));
// z0, z1 are already in the upper half-plane (positive imaginary); skip FixModuliPoint.
// Choose polarization ideals A,B so that A*B*Different(OH) is narrowly principal.
P3 := Factorization(3*OH)[1][1];   // ramified prime above 3 in Q(sqrt3) = (sqrt3)
polopts := [* <1*OH, 1*OH>, <1*OH, P3>, <P3, 1*OH>, <P3, P3>, <1*OH, P3^-1>, <P3^-1, 1*OH> *];
zlist := [* z0, z1 *];
for tag in [1..#zlist] do
    z := zlist[tag];
    LOG(Sprintf("--- moduli point %o ---", tag));
    for AB in polopts do
        try
            tau := SmallPeriodMatrix(z, AB[1], AB[2]);
            LOG(Sprintf("=== B_f PERIOD MATRIX tau (pol A=%o, B=%o) ===", AB[1], AB[2]));
            LOG(Sprintf("%o", ChangeRing(tau, ComplexField(6))));
            break;
        catch err ; end try;
    end for;
end for;
LOG("done");
exit;
