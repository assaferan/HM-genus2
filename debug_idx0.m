// debug_idx0.m -- why did the small mod-5 entries not match at weight [2,2]?
// Hypothesis: 5 is a bad prime, level-lowered at 5 => Serre weight raised.
// det constraint det = chi_5^{k0-1} = chi_5 forces k0 = 2 mod 4, so k0 in {2,6}.
// Try candidate weights on idx 0 (Q(sqrt3), p=5, serre 6) and report best match.

SetColumns(0);
d := 3; p := 5; serre := 6;
fc := [-11,2,0,5,-9,9,-3]; hc := [0,0,1];
F<w> := QuadraticField(d); OF := Integers(F);
QQ := Rationals(); Px<x> := PolynomialRing(QQ); Fp := GF(p);
fpol := &+[ fc[i]*x^(i-1) : i in [1..#fc] ];
hpol := &+[ hc[i]*x^(i-1) : i in [1..#hc] ];
C := HyperellipticCurve(fpol, hpol);
dC := Integers()!Discriminant(C);

// split-prime fingerprint of A mod 5
targ := [];
for l in PrimesInInterval(3, 400) do
    if #targ ge 20 then break; end if;
    if l eq p or (d mod l eq 0) or (dC mod l eq 0) then continue; end if;
    if not IsSquare(GF(l)!d) then continue; end if;
    fp := Factorization(l*OF);
    Append(~targ, < l, fp[1][1], fp[2][1], Integers()!(Fp!(l+1-#ChangeRing(C,GF(l)))) >);
end for;
printf "A fingerprint: %o split primes; (l,a_l mod5): %o\n", #targ, [<t[1],t[4]>: t in targ];

// norm-6 level ideal(s)
f2 := Factorization(2*OF); f3 := Factorization(3*OF);
lev := f2[1][1]*f3[1][1];  assert Norm(lev) eq 6;
printf "level norm %o\n\n", Norm(lev);

weights := [ [2,2], [2,6],[6,2], [4,6],[6,4], [6,6], [4,4], [2,4],[4,2] ];
for wt in weights do
    M := HilbertCuspForms(F, lev, wt);
    dm := Dimension(M);
    printf "weight %o: dim(new-cusp construction) %o", wt, dm;
    if dm eq 0 then printf "  (empty)\n"; continue; end if;
    D := NewformDecomposition(NewSubspace(M));
    best := 999; bestinfo := "";
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
            nb := 0;
            for pr in pairs do if red(pr[1])+red(pr[2]) ne Fq!pr[3] then nb+:=1; end if; end for;
            if nb lt best then best := nb; bestinfo := Sprintf("orbit %o deg %o resF_%o", oi, isQ select 1 else Degree(E), #Fq); end if;
        end for;
    end for;
    printf "  #orbits %o  BEST fails %o/%o  (%o)\n", #D, best, #targ, bestinfo;
end for;
exit;
