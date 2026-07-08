// match_one.m -- single (entry, weight) matcher, for the time-capped orchestrator.
// Params (Magma command line):  idx, k1, k2   e.g.  magma -b idx:=31 k1:=2 k2:=6 match_one.m
// Tests match at the small prime-to-5 Serre level, weight [k1,k2], trivial
// nebentypus, normalization t_P = Norm(P)^d * a_P(f) (d in {0,1,2}).
// Prints exactly one of:  "MATCH orbit <o> deg <D> resF_<q> d=<d>"  or  "NOMATCH".
SetColumns(0);
idx := StringToInteger(idx); k1 := StringToInteger(k1); k2 := StringToInteger(k2);
load "goal1_data.m";
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
ent := [e : e in entries | e[1] eq idx][1];
d := ent[2]; serre := ent[4]; fc := ent[5]; hc := ent[6];
F<w> := QuadraticField(d); OF := Integers(F); F5 := GF(5);
C := HyperellipticCurve(&+[fc[i]*x^(i-1):i in [1..#fc]], &+[hc[i]*x^(i-1):i in [1..#hc]]);
dC := Integers()!Discriminant(C);
J := ideal<OF|1>;
for pe in Factorization(serre) do
    P := Factorization(pe[1]*OF)[1][1]; J := J * P^(pe[2] div Valuation(Norm(P), pe[1]));
end for;
lev := J;
targ := [];
for l in PrimesInInterval(3, 300) do
    if #targ ge 14 then break; end if;
    if l eq 5 or (d mod l eq 0) or (dC mod l eq 0) or not IsSquare(GF(l)!d) then continue; end if;
    if Norm(lev) mod l eq 0 then continue; end if;
    LA := LPolynomial(ChangeRing(C, GF(l)));
    fp := Factorization(l*OF);
    Append(~targ, <l, fp[1][1], fp[2][1], -Coefficient(LA,1), Coefficient(LA,2)>);
end for;
M := HilbertCuspForms(F, lev, [k1,k2]); dm := Dimension(M);
printf "dim %o\n", dm;
if dm eq 0 then printf "NOMATCH (empty)\n"; exit; end if;
D := NewformDecomposition(NewSubspace(M));
for oi in [1..#D] do
    ef := Eigenform(D[oi]); E := HeckeEigenvalueField(D[oi]); isQ := Type(E) eq FldRat;
    lams := isQ select [* <1> *] else [* <pl[1]> : pl in Factorization(5*Integers(E)) *];
    for lam in lams do
        if isQ then Fq := F5; red := func<a|F5!(Integers()!a)>;
        else Fq,r0 := ResidueClassField(lam[1]); red := func<a|r0(Integers(E)!a)>; end if;
        if #Fq notin {5,25} then continue; end if;
        PZ := PolynomialRing(Fq); Z := PZ.1;
        for dd in [0,1,2] do
            nb := 0;
            for t in targ do
                rho := { r[1] : r in Roots(Z^2 - Fq!(F5!t[4])*Z + (Fq!(F5!t[5]) - 2*Fq!(F5!t[1]))) };
                Na := Fq!(F5!t[1]);
                if {Na^dd*red(HeckeEigenvalue(ef,t[2])), Na^dd*red(HeckeEigenvalue(ef,t[3]))} ne rho
                    then nb := 1; break; end if;
            end for;
            if nb eq 0 then printf "MATCH orbit %o deg %o resF_%o d=%o\n", oi, isQ select 1 else Degree(E), #Fq, dd; exit; end if;
        end for;
    end for;
end for;
printf "NOMATCH\n"; exit;
