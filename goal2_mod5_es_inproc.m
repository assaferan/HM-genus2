// goal2_mod5_es_inproc.m -- IN-PROCESS version of ComputePossibleModuliPoints (bypasses the
// package's ParallelPipe, which fails to spawn its L-value subprocesses here).  Replicates
// ComputeOmegaValues' assembly verbatim but computes the twisted L-values serially in-process
// via LMFDBTwistedLvalue (which works: ~10s each).  Then reconstructs B_f's period matrix and
// recognizes the genus-2 curve.
//
// Args: B (char bound, default 8), MAXN (eigenvalue norm bound), PREC (unused hint).
// Run: magma B:=8 MAXN:=10103 goal2_mod5_es_inproc.m
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
AttachSpec("../EichlerShimuraHMF/src/spec");
import "Genus2Curve.m": IgusaNumericFromTau, IgusaInvariantsInK;
if assigned B then B := StringToInteger(B); else B := 8; end if;
if assigned MAXN then MAXN := StringToInteger(MAXN); else MAXN := 10103; end if;

label := "2.2.8.1-5000.1-j";
DIR := "/scratch/home/assaferan/GitHub/HM-genus2/es_eig/";
LOGF := "goal2_mod5_es_inproc_out.txt";
PrintFile(LOGF, Sprintf("# in-process ES periods, B=%o maxn=%o", B, MAXN) : Overwrite:=true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

f := LMFDBHMFwithEigenvalues(label, DIR);
H := HeckeEigenvalueField(Parent(f));
dim := Dimension(Parent(f));
FQ := BaseField(Parent(f));
chis := QuadraticCharactersUpTo(B, FQ);
chi_signs := [HodgeSigns(chi) : chi in chis];
LOG(Sprintf("loaded; %o characters up to B=%o, dim %o", #chis, B, dim));

// --- Lvals IN-PROCESS ---
Lvals := [* *];
t0 := Cputime();
for j in [1..#chis] do
    for emb in [1..dim] do
        special := false; err := false; msg := "";
        try
            special, err := LMFDBTwistedLvalue(label, DIR, B, j, emb : maxn:=MAXN);
        catch e
            msg := "err";
        end try;
        Append(~Lvals, <j, emb, special, err, msg, 0>);
    end for;
    LOG(Sprintf("  char %o/%o done [%.1o s]", j, #chis, Cputime(t0)));
end for;

// --- assemble Omegas per sign (verbatim from ComputeOmegaValues) ---
res2 := [* *]; skipped := [* *];
desired_signs := [ [1,1], [1,-1], [-1,1], [-1,-1] ];
for s in desired_signs do
    possible_chis := Indices(chi_signs, s);
    sign_res := [* *];
    for j in possible_chis do
        chi_res := [* *];
        for emb in [1..dim] do
            entry := [x : x in Lvals | x[1] eq j and x[2] eq emb][1];
            special := entry[3]; err := entry[4]; error_message := entry[5];
            if (Type(special) eq BoolElt) or (#error_message ne 0) then continue; end if;
            prec := Precision(special);
            if Abs(special) gt 10^-(prec*0.75) then
                CC := ComplexFieldExtra(prec);
                special := CC!special;
                omega := -4*Pi(CC)^2*Sqrt(CC!2)*GaussSumOdaSimpleModuloSign(chis[j], chi_signs[j], CC)*special;
                Append(~chi_res, <omega, Abs(err)>);
            else
                Append(~skipped, entry);
            end if;
        end for;
        if #chi_res eq dim then
            Append(~sign_res, <j, <e[1] : e in chi_res>, <e[2] : e in chi_res> >);
        end if;
    end for;
    Append(~res2, sign_res);
end for;
LOG(Sprintf("omegas per sign: %o", [#sr : sr in res2]));
if 0 in [#sr : sr in res2] then
    LOG("!! a sign class has NO usable character -> cannot form moduli point (raise B or precision)");
    exit;
end if;

// --- ComputePossibleModuliPoints tail ---
Omegas_per_sign := [* [* elt[2] : elt in vps *] : vps in res2 *];
Omegas := OmegasViaCremonaTrick(H, Omegas_per_sign);
raw := PossibleModuliPoints(H, Omegas);
possible_zs := [[ FixModuliPoint(H, z) : z in zs] : zs in raw];
LOG(Sprintf("possible moduli point groups: %o", #possible_zs));

OH := Integers(H);
for gi->zs in possible_zs do
    for zi->z in zs do
        try
            tau := SmallPeriodMatrix(z, 1*OH, 1*OH);
        catch err LOG(Sprintf("  z[%o,%o]: SmallPeriodMatrix failed", gi, zi)); continue; end try;
        LOG(Sprintf("  z[%o,%o]: tau =", gi, zi));
        LOG(Sprintf("    %o", ChangeRing(tau, ComplexField(8))));
        try
            IC := IgusaNumericFromTau(tau);
            LOG(Sprintf("    IgusaNumeric = %o", [ComplexField(10)!x : x in IC]));
        catch err LOG(Sprintf("    IgusaNumericFromTau failed: %o", err`Object)); end try;
    end for;
end for;
LOG("done");
exit;
