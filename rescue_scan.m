// rescue_scan.m -- Ariel's weight/level tradeoff: can a SKIP-BIG mod-5 entry be
// matched at its SMALL prime-to-5 Serre level by RAISING the weight (Ariel found
// conductor-50000 example at level Gamma_0(4), weight [2,4])?  We scan even
// weights [k1,k2], k_i in {2,4,6}, at the Serre-conductor level, decompose into
// newforms, reduce each orbit mod lambda|5, and test the induced trace match
//   a_P(f) + a_P'(f) = a_ell(A) (mod 5)  at split primes.
//
// Run: magma rescue_scan.m
SetColumns(0);
load "goal1_data.m";
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
TARGETS := [20, 24, 34];              // smallest-serre SKIP-BIG (Q3 576, Q3 768, Q2 2304)
WEIGHTS := [[2,2],[2,4],[4,2],[2,6],[6,2],[4,4],[4,6],[6,4],[6,6]];
DIMCAP  := 9000;                       // skip a (weight) if the space is too big
NSPLIT  := 16;

LOGF := "rescue_scan_out.txt";
PrintFile(LOGF, "# rescue_scan: SKIP-BIG mod-5 at small level + raised weight" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

function IdOfNorm(OF, N)                // prime-to-p Serre ideal
    J := ideal<OF|1>;
    for pe in Factorization(N) do
        P := Factorization(pe[1]*OF)[1][1];
        J := J * P^(pe[2] div Valuation(Norm(P), pe[1]));
    end for;
    return J;
end function;

for target in TARGETS do
    ent := [e : e in entries | e[1] eq target][1];
    idx := ent[1]; d := ent[2]; p := ent[3]; serre := ent[4]; fc := ent[5]; hc := ent[6];
    F<w> := QuadraticField(d); OF := Integers(F); Fp := GF(p);
    C := HyperellipticCurve(&+[fc[i]*x^(i-1):i in [1..#fc]], &+[hc[i]*x^(i-1):i in [1..#hc]]);
    dC := Integers()!Discriminant(C);
    lev := IdOfNorm(OF, serre);
    LOG(Sprintf("\n===== idx %o  d%o p%o  Serre level norm %o =====", idx, d, p, Norm(lev)));

    // A fingerprint at split good primes
    targ := [];
    for l in PrimesInInterval(3, 400) do
        if #targ ge NSPLIT then break; end if;
        if l eq p or (d mod l eq 0) or (dC mod l eq 0) or not IsSquare(GF(l)!d) then continue; end if;
        if Norm(lev) mod l eq 0 then continue; end if;
        fp := Factorization(l*OF);
        Append(~targ, <l, fp[1][1], fp[2][1], Integers()!(Fp!(l+1-#ChangeRing(C,GF(l))))>);
    end for;

    for wt in WEIGHTS do
        t0 := Cputime();
        M := HilbertCuspForms(F, lev, wt);
        dm := Dimension(M);
        if dm eq 0 then LOG(Sprintf("  weight %o: dim 0", wt)); continue; end if;
        if dm gt DIMCAP then LOG(Sprintf("  weight %o: dim %o > cap, skip", wt, dm)); continue; end if;
        D := NewformDecomposition(NewSubspace(M));
        best := 999; bestinfo := "";
        for oi in [1..#D] do
            ef := Eigenform(D[oi]); E := HeckeEigenvalueField(D[oi]);
            isQ := Type(E) eq FldRat;
            lams := isQ select [* <1> *] else [* <pl[1]> : pl in Factorization(p*Integers(E)) *];
            pairs := [<HeckeEigenvalue(ef,t[2]), HeckeEigenvalue(ef,t[3]), t[4]> : t in targ];
            for lam in lams do
                if isQ then Fq:=Fp; red:=func<a|Fp!(Integers()!a)>;
                else Fq,r0:=ResidueClassField(lam[1]); red:=func<a|r0(Integers(E)!a)>; end if;
                nb := 0; for pr in pairs do if red(pr[1])+red(pr[2]) ne Fq!pr[3] then nb+:=1; end if; end for;
                if nb lt best then best := nb;
                   bestinfo := Sprintf("orbit %o deg %o resF_%o", oi, isQ select 1 else Degree(E), #Fq); end if;
            end for;
        end for;
        tag := (best eq 0) select "*** MATCH ***" else Sprintf("best %o/%o fail", best, #targ);
        LOG(Sprintf("  weight %o: dim %o, %o orbits -> %o (%o)  [%.1os]",
            wt, dm, #D, tag, bestinfo, Cputime(t0)));
    end for;
end for;
LOG("\nDONE.");
exit;
