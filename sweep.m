// Full sweep: for each curve match sigma to a Hilbert newform, weight [2,2].
// Level = (prime-to-2 conductor) * p2^e; sweep e from 0 up to the 2-exponent of the FULL
// conductor of A (the Serre/level-lowering upper bound), take the smallest match.
// Fingerprint t_P = -#C(F_P) mod l. Match via NewformDecomposition (dim <= DIMCAP), with a
// wrong-fingerprint discrimination control.
SetColumns(0);
load "torsion_data.m";
DIMCAP := 2600;            // skip NewformDecomposition above this (defer to kernel method)
PBOUND := 80;
outf := "sweep.out";
PrintFile(outf, "# label  d  l  base-cond  e(2-part)  sigma-level-norm  dim  orbit-dim  Hecke-deg  #primes  control" : Overwrite := true);


for row in torsion_data do
    d := row[1]; l := row[2]; lab := row[3]; fc := row[4]; hc := row[5]; N := row[6]; cc := row[7];
    K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K); Fl := GF(l);
    f := R![K| c[1]+c[2]*a : c in fc]; hh := R![K| c[1]+c[2]*a : c in hc];
    C := HyperellipticCurve(f, hh);
    p2 := Factorization(2*OK)[1][1];
    // base level (prime-to-2 conductor)
    if d eq 2 then base := ideal<OK | cc[1]+cc[2]*a>;
    else
        Cnd := Conductor(C); base := 1*OK;
        for pe in Factorization(Cnd) do if Norm(pe[1]) mod 2 ne 0 then base := base*pe[1]^pe[2]; end if; end for;
    end if;
    emax := Valuation(Conductor(C), p2);
    // fingerprint
    fp := [];
    for pp in PrimesUpTo(PBOUND) do
        if pp eq 2 or pp eq l then continue; end if;
        for tup in Factorization(pp*OK) do
            P := tup[1]; if Norm(base+P) ne 1 then continue; end if;
            k,red := ResidueClassField(P); Rk := PolynomialRing(k); ok := true; Ck := 0;
            try Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(f)], Rk![red(c):c in Coefficients(hh)]);
            catch e; ok := false; end try;
            if ok then Append(~fp, <P, Fl!(-#Ck)>); end if;
        end for;
    end for;

    matchfun := function(nf, Nlev, tw)
        eig := Eigenform(nf); Ff := HeckeEigenvalueField(nf); OF := Integers(Ff); best := 0;
        for lam in [pr[1]:pr in Factorization(l*OF)] do
            kf,m := ResidueClassField(lam); cnt := 0; okall := true;
            for pr in fp do
                P := pr[1]; t := pr[2]+tw; if Norm(Nlev+P) ne 1 then continue; end if;
                if m(OF!HeckeEigenvalue(eig,P)) ne kf!(Integers()!t) then okall := false; break; end if;
                cnt +:= 1;
            end for;
            if okall and cnt gt best then best := cnt; end if;
        end for;
        return best;
    end function;

    found := false;
    for e in [0..emax] do
        Nlev := base * p2^e;
        Mnew := NewSubspace(HilbertCuspForms(K, Nlev, [2,2])); dm := Dimension(Mnew);
        if dm eq 0 then continue; end if;
        if dm gt DIMCAP then
            line := Sprintf("%-11o d=%o l=%o base=%o e=%o levelN=%o dim=%o  DEFERRED(dim>cap)", lab,d,l,N,e,Norm(Nlev),dm);
            printf "%o\n", line; PrintFile(outf, line); found := true; break;
        end if;
        decomp := NewformDecomposition(Mnew);
        for idx in [1..#decomp] do
            mc := matchfun(decomp[idx], Nlev, 0);
            if mc ge 12 then
                ctrl := matchfun(decomp[idx], Nlev, 1);
                line := Sprintf("%-11o d=%o l=%o base=%o e=%o levelN=%o dim=%o orbit=%o Heckedeg=%o primes=%o control=%o  MATCH",
                    lab,d,l,N,e,Norm(Nlev),dm,Dimension(decomp[idx]),Degree(HeckeEigenvalueField(decomp[idx])),mc,ctrl);
                printf "%o\n", line; PrintFile(outf, line); found := true; break;
            end if;
        end for;
        if found then break; end if;
    end for;
    if not found then
        line := Sprintf("%-11o d=%o l=%o base=%o emax=%o  NO MATCH up to e_max", lab,d,l,N,emax);
        printf "%o\n", line; PrintFile(outf, line);
    end if;
end for;
PrintFile(outf, "# done"); printf "SWEEP DONE\n"; quit;
