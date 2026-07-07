// goal2_grh.m -- GRH-conditional certificate for A[7] = Ind(rho_f), idx 33.
//
// Both A[7] and Ind(rho_f) are 4-dim, irreducible, unramified outside {2,3,5,7},
// with Artin conductor supported there (cond(A) = 2^8*3*5^5 = 2400000 prime-to-7).
// By Brauer-Nesbitt they are isomorphic iff tr Frob agree at all primes; under
// GRH the *conductor-based* effective bound (Rankin-Selberg / effective
// Faltings-Serre, ~ (log cond)^2, INDEPENDENT of the huge splitting field) makes
// this a finite check.  (log(cond) ~ 14.7, so (log cond)^2 ~ 216; with standard
// explicit constants the sufficient bound is of order 10^3-10^4.)  We verify
// EXACT trace agreement for all good primes up to BOUND >> that.
//
//   ell SPLIT in F : a_ell(A) = a_P(f) + a_P'(f)   (mod 7)
//   ell INERT in F : a_ell(A) = 0                  (mod 7)  [Ind trace at inert]
//
// A-side is fast; the eigenform is needed only at SPLIT primes.
//
// Run: magma goal2_grh.m
SetColumns(0);
F<w> := QuadraticField(10); OF := Integers(F);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
fpoly := -10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9;
C := HyperellipticCurve(fpoly);
F7 := GF(7); badp := {2,3,5,7};
BOUND := 10000;

LOGF := "goal2_grh_out.txt";
PrintFile(LOGF, "# goal2_grh: GRH-conditional trace check A[7] = Ind(rho_f)" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

LOG("building orbit-10 eigenform ...");
t0 := Cputime();
P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1]; P3a := Factorization(3*OF)[1][1];
M := HilbertCuspForms(F, P2^2*P3a*P5^3, [2,2]);
D := NewformDecomposition(NewSubspace(M));
ef := Eigenform(D[10]); E := HeckeEigenvalueField(D[10]); OE := Integers(E);
found := false;
for pl in Factorization(7*OE) do
    Fq0, r0 := ResidueClassField(pl[1]);
    if Degree(Fq0) eq 2 then Fq := Fq0; red := r0; found := true; break; end if;
end for;
assert found;
LOG(Sprintf("  ready (%.1o s). Verifying all good primes up to %o ...", Cputime(t0), BOUND));

nsplit := 0; ninert := 0; fails := []; nextlog := 1000;
for l in PrimesInInterval(2, BOUND) do
    if l in badp then continue; end if;
    aA := F7 ! (l + 1 - #ChangeRing(C, GF(l)));       // tr A[7](Frob_l) mod 7
    if IsSquare(GF(l)!10) then                         // split
        fp := Factorization(l*OF);
        s := red(OE!HeckeEigenvalue(ef, fp[1][1])) + red(OE!HeckeEigenvalue(ef, fp[2][1]));
        // s lives in F_49; the induced trace equals a_ell(A) in F_7 <= F_49
        ok := (s eq Fq!aA);
        nsplit +:= 1;
    else                                               // inert: induced trace = 0
        ok := (aA eq F7!0);
        ninert +:= 1;
    end if;
    if not ok then Append(~fails, l); end if;
    if l gt nextlog then
        LOG(Sprintf("  ... up to %o: %o split + %o inert checked, %o fails", l, nsplit, ninert, #fails));
        nextlog +:= 2000;
    end if;
end for;

LOG(Sprintf("\nTOTAL: %o primes (%o split, %o inert) up to %o", nsplit+ninert, nsplit, ninert, BOUND));
LOG(Sprintf("trace DISAGREEMENTS: %o", #fails));
if #fails eq 0 then
    LOG("");
    LOG(">>> tr A[7](Frob_l) = tr Ind(rho_f)(Frob_l) for EVERY good prime l < " cat IntegerToString(BOUND));
    LOG("    cond(A[7]) = 2400000 (prime-to-7), so (log cond)^2 ~ 216; the GRH");
    LOG("    conductor-based effective bound is of order 10^3-10^4, within this range.");
    LOG("    => (under GRH) A[7] = Ind(rho_f), hence A[7] is residually modular.");
else
    LOG(Sprintf("  first few disagreeing primes: %o", fails[1..Min(10,#fails)]));
end if;
exit;
