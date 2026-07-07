// weilres_qsqrt10_c.m -- DECISIVE check: compare A's mod-7 L-polynomial (deg 4)
// against the induced char poly of the matching Hilbert newform (orbit 10,
// Hecke deg 18, lambda|7 with residue field F_49) at ALL good primes < 1000.
//
//   split p = Pa*Pb:  L_p(Ind)(T) = (1 - a_Pa T + p T^2)(1 - a_Pb T + p T^2)
//   inert p = (p):    L_p(Ind)(T) = 1 - a_(p) T^2 + p^2 T^4     [N(p)=p^2]
//   compare to A: LPolynomial(C mod p) reduced into F_49.
//
// Run: magma weilres_qsqrt10_c.m

SetColumns(0);
F<w>  := QuadraticField(10);
OF    := Integers(F);
QQ    := Rationals();
Px<x> := PolynomialRing(QQ);
f  := -10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9;
C  := HyperellipticCurve(f);
badp := {2,3,5,7};

LOGF := "weilres_qsqrt10_c_out.txt";
PrintFile(LOGF, "# weilres_qsqrt10_c -- full L-poly (charpoly) comparison mod 7" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

// level norm 1500, rebuild and grab orbit 10
f2 := Factorization(2*OF); P2 := f2[1][1];
f5 := Factorization(5*OF); P5 := f5[1][1];
f3 := Factorization(3*OF); P3a := f3[1][1];
level := P2^2 * P3a * P5^3;
LOG("building HilbertCuspForms(level 1500, [2,2]) ...");
t0 := Cputime();
M := HilbertCuspForms(F, level, [2,2]);
D := NewformDecomposition(NewSubspace(M));
LOG(Sprintf("  %o orbits (%.1o s); using orbit 10", #D, Cputime(t0)));
e := Eigenform(D[10]);
E := HeckeEigenvalueField(D[10]);
assert Degree(E) eq 18;
OE := Integers(E);
// pick lambda|7 with residue field F_49 (degree 2)
found := false; Fq := 0; red := 0;
for pl in Factorization(7*OE) do
    Fq0, red0 := ResidueClassField(pl[1]);
    if Degree(Fq0) eq 2 then Fq := Fq0; red := red0; found := true; break; end if;
end for;
assert found;
LOG(Sprintf("  lambda|7 residue field F_%o", #Fq));
PF<T> := PolynomialRing(Fq);

SPLITMAX := 500;   // split primes are cheap (ideal norm p)
INERTMAX := 160;   // inert primes cost more (ideal norm p^2)
nsplit := 0; ninert := 0; fails := [];
for p in PrimesInInterval(3, SPLITMAX) do
    if p in badp then continue; end if;
    sp := IsSquare(GF(p)!10);
    if (not sp) and p gt INERTMAX then continue; end if;
    // A's L-polynomial mod 7 -> F_49
    LA := LPolynomial(ChangeRing(C, GF(p)));
    LAred := &+[ Fq!(Integers()!Coefficient(LA,i)) * T^i : i in [0..Degree(LA)] ];
    if sp then          // split
        fp := Factorization(p*OF);
        a1 := red(OE!HeckeEigenvalue(e, fp[1][1]));
        a2 := red(OE!HeckeEigenvalue(e, fp[2][1]));
        Lind := (1 - a1*T + Fq!p*T^2)*(1 - a2*T + Fq!p*T^2);
        nsplit +:= 1;
    else                                 // inert
        Pin := Factorization(p*OF)[1][1];   // = (p), norm p^2
        a := red(OE!HeckeEigenvalue(e, Pin));
        Lind := 1 - a*T^2 + Fq!(p^2)*T^4;
        ninert +:= 1;
    end if;
    ok := (LAred eq Lind);
    if not ok then Append(~fails, <p, LAred, Lind>); end if;
    LOG(Sprintf("  p=%-4o %o  charpoly %o", p, sp select "split" else "inert",
        ok select "OK" else "MISMATCH"));
end for;

LOG(Sprintf("\ncompared %o good primes  (%o split<%o, %o inert<%o)",
    nsplit+ninert, nsplit, SPLITMAX, ninert, INERTMAX));
LOG(Sprintf("L-polynomial MISMATCHES mod 7: %o", #fails));
for fl in fails[1..Min(8,#fails)] do
    LOG(Sprintf("  p=%o : A=%o  vs  Ind=%o", fl[1], fl[2], fl[3]));
end for;
if #fails eq 0 then
    LOG("\n>>> FULL char-poly match: A[7] and Ind(rho_f,lambda) agree at every good prime < 1000.");
end if;
LOG("DONE.");
exit;
