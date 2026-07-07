// weilres_qsqrt10.m -- Ariel's mod-7 example (data idx 33).
// Curve A: y^2 = f(x), conductor 2,400,000, congruent mod 7 to a Hilbert
// newform over F = Q(sqrt(10)) induced to G_Q; Serre level norm 1500.
//
// Phase 1: (a) compute A's Frobenius traces at split primes of F (the
// discriminating data; inert primes give a_p = 0 mod 7 trivially),
// (b) build HilbertCuspForms(F, level-1500, [2,2]), decompose into newforms,
// (c) find the RATIONAL newform f with a_p(A) = a_P(f) + a_P'(f) mod 7 at all
// split test primes.
//
// Run: magma weilres_qsqrt10.m

SetColumns(0);
F<w>  := QuadraticField(10);
OF    := Integers(F);
QQ    := Rationals();
Px<x> := PolynomialRing(QQ);
F7    := GF(7);

f  := -10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9;
C  := HyperellipticCurve(f);
CondN := 2400000;            // = 2^8 * 3 * 5^5  ; bad primes {2,3,5}
badp  := {2,3,5};

LOGF := "weilres_qsqrt10_out.txt";
PrintFile(LOGF, "# weilres_qsqrt10 -- mod-7 Weil-restriction match" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

// ---- (a) A's split-prime fingerprint: <p, Pa, Pb, a_p(A) mod 7> ----
// p splits in F=Q(sqrt10) iff 10 is a QR mod p.  a_p(A) = p+1-#C(F_p).
LOG("=== A trace fingerprint at split good primes (value = a_p(A) mod 7) ===");
split_targ := [];
for p in PrimesInInterval(3, 500) do
    if p in badp then continue; end if;
    if not IsSquare(GF(p)!10) then continue; end if;   // split test
    fp := Factorization(p*OF);
    assert #fp eq 2;
    Pa := fp[1][1]; Pb := fp[2][1];
    Cp := ChangeRing(C, GF(p));
    ap := p + 1 - #Cp;
    Append(~split_targ, <p, Pa, Pb, Integers()!(F7!ap)>);
end for;
LOG(Sprintf("  %o split test primes up to 500", #split_targ));
LOG(Sprintf("  (p, a_p(A) mod 7): %o", [<t[1],t[4]> : t in split_targ[1..Min(20,#split_targ)]]));

// ---- (b) candidate level ideals of norm 1500 = 2^2 * 3 * 5^3 ----
// 2 ramifies, 5 ramifies, 3 splits in F.
f2 := Factorization(2*OF); P2 := f2[1][1];      // ramified, norm 2
f5 := Factorization(5*OF); P5 := f5[1][1];      // ramified, norm 5
f3 := Factorization(3*OF);                       // split, two primes norm 3
P3a := f3[1][1]; P3b := f3[2][1];
levels := [ P2^2 * P3a * P5^3, P2^2 * P3b * P5^3 ];
for L in levels do assert Norm(L) eq 1500; end for;
LOG(Sprintf("\n=== candidate levels of norm 1500: %o ideals ===", #levels));

// helper: does an integral eigenform e (rational Hecke field) match A mod 7?
function MatchesMod7(e, split_targ, level)
    nbad := 0; ntest := 0;
    for t in split_targ do
        p := t[1]; Pa := t[2]; Pb := t[3]; aA := t[4];
        if Norm(level) mod p eq 0 then continue; end if;   // skip bad-for-f primes
        ntest +:= 1;
        aPa := HeckeEigenvalue(e, Pa);
        aPb := HeckeEigenvalue(e, Pb);
        val := F7!(Integers()!aPa) + F7!(Integers()!aPb) - F7!aA;
        if val ne 0 then nbad +:= 1; end if;
    end for;
    return nbad, ntest;
end function;

for li -> level in levels do
    LOG(Sprintf("\n---------- level #%o, norm %o ----------", li, Norm(level)));
    t0 := Cputime();
    M := HilbertCuspForms(F, level, [2,2]);
    LOG(Sprintf("  dim(cusp) = %o  new-dim setup...", Dimension(M)));
    N := NewSubspace(M);
    D := NewformDecomposition(N);
    LOG(Sprintf("  %o newform orbits  (%.1o s)", #D, Cputime(t0)));
    for idx in [1..#D] do
        e := Eigenform(D[idx]);
        // determine Hecke field via a sample eigenvalue
        a0 := HeckeEigenvalue(e, split_targ[1][2]);
        E  := Parent(a0);
        dE := (Type(E) eq FldRat) select 1 else Degree(E);
        if dE ne 1 then
            LOG(Sprintf("  orbit %o: Hecke deg %o (non-rational, skip for E_f)", idx, dE));
            continue;
        end if;
        nbad, ntest := MatchesMod7(e, split_targ, level);
        if nbad eq 0 then
            LOG(Sprintf("  *** orbit %o: RATIONAL, MATCHES A mod 7 at %o/%o split primes ***", idx, ntest, ntest));
            // dump integral eigenvalues for E_f recovery
            evs := [];
            for t in split_targ[1..Min(12,#split_targ)] do
                Append(~evs, <Norm(t[2]), HeckeEigenvalue(e,t[2])>);
            end for;
            LOG(Sprintf("     a_P (norm, value): %o", evs));
        else
            LOG(Sprintf("  orbit %o: RATIONAL but mismatches (%o/%o fail)", idx, nbad, ntest));
        end if;
    end for;
end for;
LOG("\nDONE.");
exit;
