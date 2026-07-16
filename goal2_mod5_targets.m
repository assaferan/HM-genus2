// goal2_mod5_targets.m -- extend B_f's (s,n) = (Tr,Nm) fingerprint to split primes
// up to 200, for the CRT reconstruction (needs a larger modulus).
SetColumns(0);
d := 2; F<w> := QuadraticField(d); OF := Integers(F);
P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1];
nlev := P2^3 * P5^2;
M := HilbertCuspForms(F, nlev, [2,2]);
Dc := NewformDecomposition(NewSubspace(M));
ef := Eigenform(Dc[6]); K := HeckeEigenvalueField(Dc[6]);

LOGF := "goal2_mod5_targets_out.txt";
PrintFile(LOGF, "# split prime l : s=Tr(a_P) : n=Nm(a_P)  (B_f fingerprint)" : Overwrite := true);
for l in PrimesInInterval(3, 200) do
    if l eq 5 then continue; end if;
    if not IsSquare(GF(l)!d) then continue; end if;   // split only
    fp := Factorization(l*OF);
    aP := HeckeEigenvalue(ef, fp[1][1]);
    s := Integers()!Trace(K!aP); n := Integers()!Norm(K!aP);
    PrintFile(LOGF, Sprintf("%o %o %o", l, s, n));
    printf "l=%o s=%o n=%o\n", l, s, n;
end for;
exit;
