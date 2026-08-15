// es_eigenvalues.m -- precompute a_P(f) for the registered form 2.2.8.1-5000.1-j over
// F=Q(sqrt2) up to norm MAXN, and write them in EichlerShimuraHMF's eigenvalue-file format
// (<dir>/<label>.txt), so ComputePossibleModuliPoints can build the twisted L-functions.
//
// Args: MAXN (norm bound, default 20000).  Run: magma MAXN:=20000 es_eigenvalues.m
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
AttachSpec("../EichlerShimuraHMF/src/spec");
if assigned MAXN then MAXN := StringToInteger(MAXN); else MAXN := 20000; end if;

label := "2.2.8.1-5000.1-j";
DIR := "es_eig/";
System("mkdir -p es_eig");

LOGF := "es_eigenvalues_out.txt";
PrintFile(LOGF, Sprintf("# eigenvalues for %o up to norm %o", label, MAXN) : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

f := LMFDBHMF(label);
Mf := Parent(f);
F := BaseField(Mf); OF := Integers(F);
LOG(Sprintf("form loaded: dim %o, Hecke field %o", Dimension(Mf), DefiningPolynomial(HeckeEigenvalueField(Mf))));

// initialize + RESUME: load any already-computed eigenvalues so we only fill the gaps
if not assigned f`hecke_eigenvalues then f`hecke_eigenvalues := AssociativeArray(); end if;
efile := DIR cat label cat ".txt";
if OpenTest(efile, "r") then
    LoadEigenvalues(~f, efile);
    LOG(Sprintf("resumed: %o eigenvalues already present", #Keys(f`hecke_eigenvalues)));
end if;

primes := [ P : P in PrimesUpTo(MAXN, F) | not IsDefined(f`hecke_eigenvalues, P) ];
LOG(Sprintf("#primes up to norm %o still to compute = %o", MAXN, #primes));
t0 := Cputime(); cnt := 0; failed := [];
for P in primes do
    try
        f`hecke_eigenvalues[P] := HeckeEigenvalue(f, P);
    catch err
        Append(~failed, P);
        LOG(Sprintf("  !! HeckeEigenvalue FAILED at prime norm %o (%o): %o", Norm(P), LMFDBLabel(P), err`Object));
    end try;
    cnt +:= 1;
    if cnt mod 250 eq 0 then
        LOG(Sprintf("  %o/%o primes, %.1o s (%.3o s/prime)", cnt, #primes, Cputime(t0), Cputime(t0)/cnt));
        _ := WriteEigenvalues(f, DIR cat label cat ".txt" : Overwrite := true);  // checkpoint
    end if;
end for;
LOG(Sprintf("computed %o eigenvalues in %.1o s (%o failed: %o)", cnt-#failed, Cputime(t0), #failed, [Norm(P) : P in failed]));

n := WriteEigenvalues(f, DIR cat label cat ".txt" : Overwrite := true);
LOG(Sprintf("wrote %o eigenvalues to %o%o.txt", n, DIR, label));
// report the effective norm bound the package will see
bnd := NormBoundOnComputedEigenvalues(f);
LOG(Sprintf("NormBoundOnComputedEigenvalues = %o", bnd));
exit;
