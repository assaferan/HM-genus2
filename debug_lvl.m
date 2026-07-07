// debug_lvl.m -- confirm the mod-5 recipe: level = SerreIdeal * P5^e, weight [2,2].
SetColumns(0);
load "goal1_data.m";
QQ := Rationals(); Px<x> := PolynomialRing(QQ);

function IdOfNorm(OF, N)  // one ideal of norm N (prime-to-p part, squarefree-ish here)
    J := ideal<OF|1>;
    for pe in Factorization(N) do
        Ps := [pl[1] : pl in Factorization(pe[1]*OF)];
        // greedily use first prime; exponent so that norm matches
        P := Ps[1]; fdeg := Valuation(Norm(P), pe[1]);
        assert pe[2] mod fdeg eq 0;
        J := J * P^(pe[2] div fdeg);
    end for;
    return J;
end function;

for target in [5, 8] do
    ent := [e : e in entries | e[1] eq target][1];
    idx := ent[1]; d := ent[2]; p := ent[3]; serre := ent[4]; fc := ent[5]; hc := ent[6];
    F<w> := QuadraticField(d); OF := Integers(F); Fp := GF(p);
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
    base := IdOfNorm(OF, serre);
    P5 := Factorization(p*OF)[1][1];
    printf "\n===== idx %o  d%o p%o serre%o  (NormP5=%o) =====\n", idx, d, p, serre, Norm(P5);
    for e in [1,2] do
        lev := base * P5^e;
        M := HilbertCuspForms(F, lev, [2,2]); dm := Dimension(M);
        printf "  P5^%o: level norm %o, dim %o", e, Norm(lev), dm;
        if dm eq 0 then printf "  empty\n"; continue; end if;
        D := NewformDecomposition(NewSubspace(M));
        nmatch := 0; bestfail := 999;
        for oi in [1..#D] do
            ef := Eigenform(D[oi]); E := HeckeEigenvalueField(D[oi]);
            isQ := Type(E) eq FldRat;
            lams := isQ select [* <1> *] else [* <pl[1]> : pl in Factorization(p*Integers(E)) *];
            pairs := [<HeckeEigenvalue(ef,t[2]),HeckeEigenvalue(ef,t[3]),t[4]> : t in targ | Norm(lev) mod t[1] ne 0];
            for lam in lams do
                if isQ then Fq:=Fp; red:=func<a|Fp!(Integers()!a)>;
                else Fq,r0:=ResidueClassField(lam[1]); red:=func<a|r0(Integers(E)!a)>; end if;
                nb := 0; for pr in pairs do if red(pr[1])+red(pr[2]) ne Fq!pr[3] then nb+:=1; end if; end for;
                if nb eq 0 then nmatch +:= 1; end if;
                if nb lt bestfail then bestfail := nb; end if;
            end for;
        end for;
        printf "  #orbits %o  matches %o  bestfail %o/%o\n", #D, nmatch, bestfail, #targ;
    end for;
end for;
exit;
