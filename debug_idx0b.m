// debug_idx0b.m -- det forces weight [2,2]; Serre cond 6 is prime-to-5.
// Since 5 is a bad prime (inert in Q(sqrt3), norm 25), restore the 5-part in
// the LEVEL and test weight [2,2] at norm 6, 6*25, 6*25^2.

SetColumns(0);
d := 3; p := 5;
fc := [-11,2,0,5,-9,9,-3]; hc := [0,0,1];
F<w> := QuadraticField(d); OF := Integers(F);
QQ := Rationals(); Px<x> := PolynomialRing(QQ); Fp := GF(p);
C := HyperellipticCurve(&+[fc[i]*x^(i-1):i in [1..#fc]], &+[hc[i]*x^(i-1):i in [1..#hc]]);
dC := Integers()!Discriminant(C);

targ := [];
for l in PrimesInInterval(3, 500) do
    if #targ ge 24 then break; end if;
    if l eq p or (d mod l eq 0) or (dC mod l eq 0) then continue; end if;
    if not IsSquare(GF(l)!d) then continue; end if;
    fp := Factorization(l*OF);
    Append(~targ, <l, fp[1][1], fp[2][1], Integers()!(Fp!(l+1-#ChangeRing(C,GF(l))))>);
end for;

P2 := Factorization(2*OF)[1][1];
P3 := Factorization(3*OF)[1][1];
P5 := Factorization(5*OF)[1][1];   // 5 inert, norm 25
base := P2*P3;
printf "Norm(P5)=%o\n", Norm(P5);

for lev in [ base, base*P5, base*P5^2 ] do
    M := HilbertCuspForms(F, lev, [2,2]); dm := Dimension(M);
    printf "\nlevel norm %o (=%o): dim %o\n", Norm(lev), Factorization(Norm(lev)), dm;
    if dm eq 0 then printf "  empty\n"; continue; end if;
    D := NewformDecomposition(NewSubspace(M));
    printf "  %o orbits\n", #D;
    for oi in [1..#D] do
        e := Eigenform(D[oi]); E := HeckeEigenvalueField(D[oi]);
        isQ := Type(E) eq FldRat;
        lams := isQ select [* <1> *] else [* <pl[1]> : pl in Factorization(p*Integers(E)) *];
        pairs := [];
        for t in targ do
            if Norm(lev) mod t[1] eq 0 then continue; end if;
            Append(~pairs, <HeckeEigenvalue(e,t[2]), HeckeEigenvalue(e,t[3]), t[4]>);
        end for;
        for lam in lams do
            if isQ then Fq := Fp; red := func<a|Fp!(Integers()!a)>;
            else Fq,r0 := ResidueClassField(lam[1]); red := func<a|r0(Integers(E)!a)>; end if;
            nb := 0; for pr in pairs do if red(pr[1])+red(pr[2]) ne Fq!pr[3] then nb+:=1; end if; end for;
            tag := (nb eq 0) select "*** MATCH ***" else Sprintf("%o/%o fail", nb, #pairs);
            printf "    orbit %o deg %o resF_%o: %o\n", oi, isQ select 1 else Degree(E), #Fq, tag;
        end for;
    end for;
end for;
exit;
