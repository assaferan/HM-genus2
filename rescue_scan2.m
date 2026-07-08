// rescue_scan2.m -- CORRECTED weight-raising rescue.  Small (prime-to-5) Serre
// level, raised weight, TRIVIAL nebentypus, with the normalization
//   t_P = Norm(P)^d * a_P(f) mod lambda      (d = weight-dependent twist; try 0,1,2)
// matched against the two sigma-traces of A at split primes.  (Fixes the first
// scan, which compared raw a_P(f); Ariel's [2,4] match uses d=1.)
//
// Run: magma rescue_scan2.m
SetColumns(0);
load "goal1_data.m";
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
TARGETS := [8, 0, 20, 24, 34];       // 8=in-data validation; 0=known niveau-2; 20/24/34=SKIP-BIG
WEIGHTS := [[2,2],[2,4],[4,2],[2,6],[6,2],[4,4]];
DVALS   := [0,1,2];
DIMCAP  := 3000;
NSPLIT  := 14;

LOGF := "rescue_scan2_out.txt";
PrintFile(LOGF, "# rescue_scan2: small level + raised weight + Norm(P)^d normalization" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

function IdOfNorm(OF, N)
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
    F<w> := QuadraticField(d); OF := Integers(F); F5 := GF(5); Fp := F5;
    C := HyperellipticCurve(&+[fc[i]*x^(i-1):i in [1..#fc]], &+[hc[i]*x^(i-1):i in [1..#hc]]);
    dC := Integers()!Discriminant(C);
    lev := IdOfNorm(OF, serre);
    LOG(Sprintf("\n===== idx %o  d%o  Serre level norm %o =====", idx, d, Norm(lev)));

    // A: split primes -> store <l, Pa, Pb, a1, b1> (A-side integers); the two
    // sigma-traces (rho) are computed later INSIDE the eigenform's residue field.
    targ := [];
    for l in PrimesInInterval(3, 300) do
        if #targ ge NSPLIT then break; end if;
        if l eq p or (d mod l eq 0) or (dC mod l eq 0) or not IsSquare(GF(l)!d) then continue; end if;
        if Norm(lev) mod l eq 0 then continue; end if;
        LA := LPolynomial(ChangeRing(C, GF(l)));
        fp := Factorization(l*OF);
        Append(~targ, <l, fp[1][1], fp[2][1], -Coefficient(LA,1), Coefficient(LA,2)>);
    end for;

    for wt in WEIGHTS do
        M := HilbertCuspForms(F, lev, wt); dm := Dimension(M);
        if dm eq 0 then continue; end if;
        if dm gt DIMCAP then LOG(Sprintf("  weight %o: dim %o > cap, skip", wt, dm)); continue; end if;
        D := NewformDecomposition(NewSubspace(M));
        best := 999; binfo := "";
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
                    nb := 0;
                    for t in targ do
                        // rho = the two sigma-traces of A, as roots over Fq
                        rho := { r[1] : r in Roots(Z^2 - Fq!(F5!t[4])*Z + (Fq!(F5!t[5]) - 2*Fq!(F5!t[1]))) };
                        Na := Fq!(F5!(t[1]));    // Norm(Pa)=l
                        va := (Na^dd) * red(HeckeEigenvalue(ef, t[2]));
                        vb := (Na^dd) * red(HeckeEigenvalue(ef, t[3]));
                        if {va,vb} ne rho then nb +:= 1; end if;
                    end for;
                    if nb lt best then best := nb;
                       binfo := Sprintf("orbit %o deg %o resF_%o d=%o", oi, isQ select 1 else Degree(E), #Fq, dd); end if;
                end for;
            end for;
        end for;
        tag := (best eq 0) select "*** MATCH ***" else Sprintf("best %o/%o", best, #targ);
        LOG(Sprintf("  weight %o: dim %o, %o orbits -> %o (%o)", wt, dm, #D, tag, binfo));
    end for;
end for;
LOG("\nDONE.");
exit;
