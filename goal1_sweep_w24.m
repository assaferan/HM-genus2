// goal1_sweep_w24.m -- Goal 1, mod-5, DONE RIGHT: match at the SMALL prime-to-5
// Serre level with RAISED weight [2,4] (Ariel's approach), trivial nebentypus,
// normalization t_P = Norm(P)^d * a_P(f) mod lambda (d in {0,1,2}; d=1 for [2,4]).
// This avoids the ×625 blow-up of the weight-[2,2] approach, so ALL 70 mod-5
// entries (incl. the former 50 "SKIP-BIG") are feasible.
//
// Resume-aware; incremental.  Run: magma goal1_sweep_w24.m
SetColumns(0);
load "goal1_data.m";
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
WEIGHTS := [[2,4],[4,2]];
DVALS   := [1,0,2];
DIMCAP  := 22000;
NSPLIT  := 16;
RESF    := "goal1_w24_results.txt";

done := {};
if OpenTest(RESF, "r") then
    fh := Open(RESF, "r");
    while true do s := Gets(fh); if IsEof(s) then break; end if;
        if #s ge 5 and s[1..5] eq "ENTRY" and Position(s,"SKIP-BIG") eq 0 then
            done := done join {StringToInteger(Split(s," \t")[2])}; end if;
    end while; delete fh;
else
    PrintFile(RESF, "# Goal-1 mod-5 at small level + weight [2,4]: ENTRY idx d serre -> status");
end if;
procedure REC(s) printf "%o\n", s; PrintFile(RESF, s); end procedure;

function IdOfNorm(OF, N)
    J := ideal<OF|1>;
    for pe in Factorization(N) do
        P := Factorization(pe[1]*OF)[1][1];
        J := J * P^(pe[2] div Valuation(Norm(P), pe[1]));
    end for;
    return J;
end function;

for ent in entries do
    idx := ent[1]; d := ent[2]; p := ent[3]; serre := ent[4]; fc := ent[5]; hc := ent[6];
    if p ne 5 then continue; end if;           // mod-5 only
    if idx in done then continue; end if;
    F<w> := QuadraticField(d); OF := Integers(F); F5 := GF(5);
    C := HyperellipticCurve(&+[fc[i]*x^(i-1):i in [1..#fc]], &+[hc[i]*x^(i-1):i in [1..#hc]]);
    dC := Integers()!Discriminant(C);
    lev := IdOfNorm(OF, serre);
    hdr := Sprintf("ENTRY %o d%o p5 serre%o (Serre level norm %o)", idx, d, serre, Norm(lev));

    // A-side split-prime data
    targ := [];
    for l in PrimesInInterval(3, 400) do
        if #targ ge NSPLIT then break; end if;
        if l eq 5 or (d mod l eq 0) or (dC mod l eq 0) or not IsSquare(GF(l)!d) then continue; end if;
        if Norm(lev) mod l eq 0 then continue; end if;
        LA := LPolynomial(ChangeRing(C, GF(l)));
        fp := Factorization(l*OF);
        Append(~targ, <l, fp[1][1], fp[2][1], -Coefficient(LA,1), Coefficient(LA,2)>);
    end for;

    matched := false; t0 := Cputime(); skipbig := false;
    for wt in WEIGHTS do
        M := HilbertCuspForms(F, lev, wt); dm := Dimension(M);
        if dm eq 0 then continue; end if;
        if dm gt DIMCAP then skipbig := true; continue; end if;
        D := NewformDecomposition(NewSubspace(M));
        for oi in [1..#D] do
            ef := Eigenform(D[oi]); E := HeckeEigenvalueField(D[oi]);
            isQ := Type(E) eq FldRat;
            lams := isQ select [* <1> *] else [* <pl[1]> : pl in Factorization(5*Integers(E)) *];
            for lam in lams do
                if isQ then Fq := F5; red := func<a|F5!(Integers()!a)>;
                else Fq,r0 := ResidueClassField(lam[1]); red := func<a|r0(Integers(E)!a)>; end if;
                if #Fq notin {5,25} then continue; end if;
                PZ := PolynomialRing(Fq); Z := PZ.1;
                for dd in DVALS do
                    nb := 0; nt := 0;
                    for t in targ do
                        nt +:= 1;
                        rho := { r[1] : r in Roots(Z^2 - Fq!(F5!t[4])*Z + (Fq!(F5!t[5]) - 2*Fq!(F5!t[1]))) };
                        Na := Fq!(F5!(t[1]));
                        va := (Na^dd) * red(HeckeEigenvalue(ef, t[2]));
                        vb := (Na^dd) * red(HeckeEigenvalue(ef, t[3]));
                        if {va,vb} ne rho then nb +:= 1; end if;
                    end for;
                    if nb eq 0 and nt ge 10 then
                        REC(Sprintf("%o -> MATCH weight %o orbit %o HeckeDeg %o resF_%o d=%o tests %o [%.1os]",
                            hdr, wt, oi, isQ select 1 else Degree(E), #Fq, dd, nt, Cputime(t0)));
                        matched := true; break;
                    end if;
                end for;
                if matched then break; end if;
            end for;
            if matched then break; end if;
        end for;
        if matched then break; end if;
    end for;
    if not matched then
        if skipbig then REC(Sprintf("%o -> SKIP-BIG (dim > %o at [2,4])", hdr, DIMCAP));
        else REC(Sprintf("%o -> NO MATCH [%.1os]", hdr, Cputime(t0))); end if;
    end if;
end for;
REC("# pass done");
exit;
