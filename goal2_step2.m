// goal2_step2.m -- Goal 2, Step 2 (GL2/F approach) for idx 33.
// Verify sigma_A = rho_f as 2-dim mod-7 reps of G_F = G_{Q(sqrt10)} by matching
// TRACES at primes of F (the empirical input to a Faltings-Serre certificate):
//   * split l = P*P' : {t1,t2} from A (roots of Z^2 - a1 Z + (a2-2l)) must equal
//     {a_P(f), a_P'(f)} mod lambda, as multisets in F_49.
//   * inert l = (l)   : tr sigma_A(Frob_(l)) = -a2(A,l) mod 7 must equal
//     a_(l)(f) mod lambda.  (Unambiguous single-valued test; the strong one.)
// det on both sides = Norm(P) mod 7 (weight 2), automatic.
//
// Run: magma goal2_step2.m

SetColumns(0);
F<w> := QuadraticField(10); OF := Integers(F);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
f := -10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9;
C := HyperellipticCurve(f);
F7 := GF(7); badp := {2,3,5,7};

LOGF := "goal2_step2_out.txt";
PrintFile(LOGF, "# goal2_step2: sigma_A = rho_f over F=Q(sqrt10)" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

// ---- rebuild orbit-10 eigenform, lambda|7 with residue field F_49 ----
LOG("building HilbertCuspForms(1500,[2,2]) and locating orbit 10 ...");
t0 := Cputime();
P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1];
P3a := Factorization(3*OF)[1][1];
level := P2^2 * P3a * P5^3;
M := HilbertCuspForms(F, level, [2,2]);
D := NewformDecomposition(NewSubspace(M));
e := Eigenform(D[10]); E := HeckeEigenvalueField(D[10]); OE := Integers(E);
assert Degree(E) eq 18;
found := false;
for pl in Factorization(7*OE) do
    Fq0, r0 := ResidueClassField(pl[1]);
    if Degree(Fq0) eq 2 then Fq := Fq0; red := r0; found := true; break; end if;
end for;
assert found;
LOG(Sprintf("  orbit 10 ready, residue field F_%o  (%.1o s)", #Fq, Cputime(t0)));
PZ<Z> := PolynomialRing(Fq);

// ---- SPLIT primes of F: multiset match ----
nsplit := 0; splitfail := 0;
for l in PrimesInInterval(3, 400) do
    if l in badp or not IsSquare(GF(l)!10) then continue; end if;
    LA := LPolynomial(ChangeRing(C, GF(l)));
    a1 := -Coefficient(LA,1); a2 := Coefficient(LA,2);
    // sigma_A traces {t1,t2} = roots of Z^2 - a1 Z + (a2 - 2l) over Fq
    tset := { rt[1] : rt in Roots(Z^2 - Fq!(F7!a1)*Z + (Fq!(F7!a2) - 2*Fq!(F7!l))) };
    fp := Factorization(l*OF);
    aset := { red(OE!HeckeEigenvalue(e, fp[1][1])), red(OE!HeckeEigenvalue(e, fp[2][1])) };
    if tset ne aset then splitfail +:= 1; end if;
    nsplit +:= 1;
end for;
LOG(Sprintf("SPLIT primes of F: %o tested, %o multiset mismatches", nsplit, splitfail));

// ---- INERT primes of F: direct single-valued match (slow; small l only) ----
ninert := 0; inertfail := 0; details := [];
for l in PrimesInInterval(3, 45) do
    if l in badp or IsSquare(GF(l)!10) then continue; end if;   // want inert
    LA := LPolynomial(ChangeRing(C, GF(l)));
    trA := Fq!(F7!(-Coefficient(LA,2)));                          // -a2 mod 7
    Pin := Factorization(l*OF)[1][1];                             // = (l), norm l^2
    trf := red(OE!HeckeEigenvalue(e, Pin));
    ok := (trA eq trf);
    if not ok then inertfail +:= 1; end if;
    Append(~details, <l, trA, trf, ok>);
    ninert +:= 1;
end for;
LOG(Sprintf("INERT primes of F (l<=45): %o tested, %o mismatches", ninert, inertfail));
for dd in details do
    LOG(Sprintf("   l=%o (norm %o): tr sigma_A=%o  a_(l)(f)=%o  %o",
        dd[1], dd[1]^2, dd[2], dd[3], dd[4] select "OK" else "FAIL"));
end for;

LOG("");
if splitfail eq 0 and inertfail eq 0 then
    LOG(">>> sigma_A and rho_f have identical traces at ALL tested primes of F");
    LOG("    (split multisets + inert single values).  This is the trace-agreement");
    LOG("    input for the Faltings-Serre certificate over F; combined with big image");
    LOG("    (Step 1) it drives sigma_A = rho_f, hence A[7] = Ind(rho_f).");
end if;
exit;
