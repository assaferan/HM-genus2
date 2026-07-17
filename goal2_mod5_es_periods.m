// goal2_mod5_es_periods.m -- compute B_f's period matrix via EichlerShimuraHMF (Oda periods
// from twisted L-values) and reconstruct the genus-2 curve using this repo's own tau->curve.
//
// ComputePossibleModuliPoints -> SmallPeriodMatrix (2x2 tau) -> IgusaNumericFromTau ->
// (try) recognize Igusa invariants in Q(sqrt2).  B_f is dim 2 => no Schottky needed.
//
// Args: B (char bound, default 15), MAXN (eigenvalue norm bound; default false=auto),
//       CORES (default 8).  Run: magma B:=15 MAXN:=2999 CORES:=8 goal2_mod5_es_periods.m
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
AttachSpec("../EichlerShimuraHMF/src/spec");
import "Genus2Curve.m": IgusaNumericFromTau, IgusaInvariantsInK;
if assigned B then B := StringToInteger(B); else B := 15; end if;
if assigned CORES then CORES := StringToInteger(CORES); else CORES := 8; end if;
if assigned MAXN then MAXN := StringToInteger(MAXN); else MAXN := false; end if;

label := "2.2.8.1-5000.1-j";
DIR := "/scratch/home/assaferan/GitHub/HM-genus2/es_eig/";  // absolute: subprocesses may run elsewhere
LOGF := "goal2_mod5_es_periods_out.txt";
PrintFile(LOGF, Sprintf("# ES periods for %o, B=%o maxn=%o", label, B, MAXN) : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

f := LMFDBHMFwithEigenvalues(label, DIR);
OH := Integers(HeckeEigenvalueField(Parent(f)));
FQ := BaseField(Parent(f));                 // Q(sqrt2)
LOG(Sprintf("form + eigenvalues loaded; effective norm bound %o", NormBoundOnComputedEigenvalues(f)));

t0 := Cputime();
possible_zs := ComputePossibleModuliPoints(CORES, label, DIR, B : maxn := MAXN);
LOG(Sprintf("ComputePossibleModuliPoints -> %o groups [%.1o s]", #possible_zs, Cputime(t0)));

nseen := 0;
for gi->zs in possible_zs do
    for zi->z in zs do
        nseen +:= 1;
        try
            tau := SmallPeriodMatrix(z, 1*OH, 1*OH);
        catch err
            LOG(Sprintf("  z[%o,%o]: SmallPeriodMatrix failed: %o", gi, zi, err`Object)); continue;
        end try;
        // diagnostic: Schottky (should be tiny for a genuine period matrix), then Igusa
        schok := true; sch := 0;
        try sch := Abs(EvaluateSchottkyModularForm(tau)); catch err schok := false; end try;
        LOG(Sprintf("  z[%o,%o]: tau ok; |Schottky|=%o", gi, zi, schok select RealField(6)!sch else "n/a"));
        try
            IC := IgusaNumericFromTau(tau);           // numeric absolute Igusa invariants
            LOG(Sprintf("    IgusaNumeric = %o", [ComplexField(8)!x : x in IC]));
        catch err
            LOG(Sprintf("    IgusaNumericFromTau failed: %o", err`Object)); continue;
        end try;
        // try to recognize the invariants in Q(sqrt2) (both embeddings)
        for e in InfinitePlaces(FQ) do
            emb := hom< FQ -> Universe(IC) | Evaluate(FQ.1, e : Precision := Precision(Universe(IC))) >;
            try
                IK := IgusaInvariantsInK(tau, FQ, emb);
                LOG(Sprintf("    *** recognized Igusa in Q(sqrt2) [emb %o]: %o ***", e, IK));
            catch err
                LOG(Sprintf("    recognize in Q(sqrt2) [emb %o] failed (precision?): %o", e, err`Object));
            end try;
        end for;
    end for;
end for;
LOG(Sprintf("done: %o moduli points examined", nseen));
exit;
