// weilres_qsqrt10_b.m -- test ALL newform orbits at level norm 1500 for a
// mod-7 congruence with A[7] (induced).  f need NOT be rational: we reduce each
// orbit's eigenvalues mod every prime lambda | 7 of its Hecke field and test
//   a_p(A) = a_Pa(f) + a_Pb(f)   in the residue field, at split primes p.
//
// Run: magma weilres_qsqrt10_b.m

SetColumns(0);
F<w>  := QuadraticField(10);
OF    := Integers(F);
QQ    := Rationals();
Px<x> := PolynomialRing(QQ);

f  := -10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9;
C  := HyperellipticCurve(f);
badp := {2,3,5};

LOGF := "weilres_qsqrt10_b_out.txt";
PrintFile(LOGF, "# weilres_qsqrt10_b -- mod-7 reduction of all orbits" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

// A's split-prime fingerprint: <p, Pa, Pb, a_p(A) in Z>
split_targ := [];
for p in PrimesInInterval(3, 500) do
    if p in badp then continue; end if;
    if not IsSquare(GF(p)!10) then continue; end if;
    fp := Factorization(p*OF); assert #fp eq 2;
    Cp := ChangeRing(C, GF(p));
    Append(~split_targ, <p, fp[1][1], fp[2][1], p+1-#Cp>);
end for;
LOG(Sprintf("%o split test primes up to 500", #split_targ));

// level #1 of norm 1500
f2 := Factorization(2*OF); P2 := f2[1][1];
f5 := Factorization(5*OF); P5 := f5[1][1];
f3 := Factorization(3*OF); P3a := f3[1][1];
level := P2^2 * P3a * P5^3; assert Norm(level) eq 1500;
LOG(Sprintf("level norm %o", Norm(level)));

t0 := Cputime();
M := HilbertCuspForms(F, level, [2,2]);
D := NewformDecomposition(NewSubspace(M));
LOG(Sprintf("%o orbits (%.1o s)\n", #D, Cputime(t0)));

for idx in [1..#D] do
    e  := Eigenform(D[idx]);
    a0 := HeckeEigenvalue(e, split_targ[1][2]);
    E  := Parent(a0);
    if Type(E) eq FldRat then
        lams := [* <1, "Q", func<a| GF(7)!(Integers()!a)>, GF(7)> *];
    else
        OE := Integers(E);
        lams := [* *];
        for pl in Factorization(7*OE) do
            lam := pl[1]; Fq, red := ResidueClassField(lam);
            Append(~lams, <Degree(Fq), Sprintf("deg%o", Degree(Fq)),
                           func<a| red(OE!a)>, Fq>);
        end for;
    end if;
    dEs := (Type(E) eq FldRat) select 1 else Degree(E);
    for lam in lams do
        redf := lam[3]; Fq := lam[4];
        nbad := 0; ntest := 0;
        for t in split_targ do
            p := t[1];
            if Norm(level) mod p eq 0 then continue; end if;
            ntest +:= 1;
            va := redf(HeckeEigenvalue(e, t[2]));
            vb := redf(HeckeEigenvalue(e, t[3]));
            if va + vb ne Fq!(t[4]) then nbad +:= 1; end if;
        end for;
        tag := (nbad eq 0) select "*** MATCH ***" else Sprintf("%o/%o fail", nbad, ntest);
        LOG(Sprintf("orbit %2o (Hecke deg %2o), lambda|7 res %o: %o",
            idx, dEs, lam[2], tag));
    end for;
end for;
LOG("\nDONE.");
exit;
