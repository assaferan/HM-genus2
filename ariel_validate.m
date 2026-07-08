// ariel_validate.m -- reproduce Ariel's mod-5 match (mod-5-match.txt) to fix the
// normalization: at raised weight [2,4] the match uses t_P = Norm(P)*a_P(f) mod
// lambda, NOT raw a_P(f), with TRIVIAL nebentypus.
SetColumns(0);
F<w> := QuadraticField(2); OF := Integers(F);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
// C: y^2 + x^3 y = -2x^6 -5x^5 -5x^4 +5x^2 +4x
C := HyperellipticCurve(-2*x^6-5*x^5-5*x^4+5*x^2+4*x, x^3);
P2 := Factorization(2*OF)[1][1];               // ramified, norm 2
level := P2^4;  assert Norm(level) eq 16;
t0 := Cputime();
M := HilbertCuspForms(F, level, [2,4]);
printf "space built, dim %o (%.1os)\n", Dimension(M), Cputime(t0); t0:=Cputime();
D := NewformDecomposition(NewSubspace(M));
printf "decomposed, %o orbits (%.1os)\n", #D, Cputime(t0);
ef := Eigenform(D[1]); E := HeckeEigenvalueField(D[1]);
printf "Hecke field degree %o; #orbits %o\n", (Type(E) eq FldRat) select 1 else Degree(E), #D;
OE := Integers(E);
lam := Factorization(5*OE)[1][1]; Fq, red := ResidueClassField(lam);
printf "residue field F_%o\n\n", #Fq;
PZ<Z> := PolynomialRing(Fq);

printf "  p  type   (a_p,b_p)    rho(from A)        t=Norm*a_P(f)      ok\n";
nok := 0; ntot := 0;
for p in PrimesInInterval(3, 113) do
    if p in {2,5} then continue; end if;
    LA := LPolynomial(ChangeRing(C, GF(p)));
    ap := Integers()!(-Coefficient(LA,1)); bp := Integers()!Coefficient(LA,2);
    split := IsSquare(GF(p)!2);
    if split then
        // rho traces from A: roots of Z^2 - ap Z + (bp - 2p)  (the two sigma-traces)
        rho := { r[1] : r in Roots(Z^2 - Fq!(GF(5)!ap)*Z + (Fq!(GF(5)!bp) - 2*Fq!(GF(5)!p))) };
        fp := Factorization(p*OF);
        tset := { Fq!(GF(5)!p) * red(OE!HeckeEigenvalue(ef, fp[1][1])),
                  Fq!(GF(5)!p) * red(OE!HeckeEigenvalue(ef, fp[2][1])) };
        ok := (rho eq tset); typ := "split";
    else
        rho := { Fq!(GF(5)!(-bp)) };            // inert: sigma-trace = -b_p
        Pin := Factorization(p*OF)[1][1];       // norm p^2
        tset := { Fq!(GF(5)!(p^2)) * red(OE!HeckeEigenvalue(ef, Pin)) };
        ok := (rho eq tset); typ := "inert";
    end if;
    ntot +:= 1; if ok then nok +:= 1; end if;
    printf "%3o  %5o  (%o,%o)  %-18o %-18o %o\n", p, typ, ap, bp,
        Sprint([x : x in rho]), Sprint([x : x in tset]), ok select "yes" else "NO";
end for;
printf "\nMATCH: %o / %o primes\n", nok, ntot;
if nok eq ntot then printf ">>> Ariel's [2,4]/level-16 match REPRODUCED with t_P = Norm(P)*a_P(f), trivial nebentypus.\n"; end if;
exit;
