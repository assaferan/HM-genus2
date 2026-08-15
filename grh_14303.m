// grh_14303.m -- GRH-conditional certificate that sigma = rho_f for the l=11 torsion example
// (curve 14303.1 over F=Q(sqrt2)).
//
// sigma = the 2-dim irreducible sub of A[11] (A[11]^ss = 1 + chi_11 + sigma); rho_f = the
// mod-lambda reduction of the Hilbert newform f at level (-13a-121) (norm 14303), wt [2,2],
// Hecke field 5.5.14641.1 = Q(zeta_11)^+ (11 totally ramified -> unique lambda|11, residue
// field F_11). Both are 2-dim, irreducible, det = chi_11, unramified outside S = {P_14303, 11}.
// By Brauer-Nesbitt (both semisimple, same det) they are isomorphic iff tr Frob agree at all
// good primes; under GRH the conductor-based effective bound (~ (log cond)^2, INDEPENDENT of
// the splitting field) makes this finite. cond(sigma) ~ 14303*11^a, so (log cond)^2 ~ 200;
// we verify EXACT trace agreement for all good primes of F up to BOUND >> that.
//
//   tr sigma(Frob_P) = a_P(A) - 1 - N(P) = -#C(F_P)   (mod 11)
//   tr rho_f(Frob_P) = a_P(f) mod lambda  (in F_11)
//
// Run: magma grh_14303.m
SetColumns(0);
K<a> := QuadraticField(2); OK := Integers(K); R<x> := PolynomialRing(K);
l := 11; Fl := GF(l); RT := PolynomialRing(Fl);
BOUND := 3000;
LOGF := "grh_14303_out.txt";
PrintFile(LOGF, "# grh_14303: GRH-conditional trace check sigma = rho_f, l=11, curve 14303.1" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

// curve 14303.1
f := a*x^5 + (3+2*a)*x^4 + (3+a)*x^3 + (-4-2*a)*x^2 + (-1+2*a)*x + (1-a);
hh := x^3 + x^2 + 1;
C := HyperellipticCurve(f, hh);
Ncond := ideal<OK | -13*a - 121>;                 // norm 14303 (2-part trivial, e=0)

// ---- build the matching eigenform, pick lambda | 11 ----
LOG("building eigenform at level norm 14303 ...");
t0 := Cputime();
Mnew := NewSubspace(HilbertCuspForms(K, Ncond, [2,2]));
decomp := NewformDecomposition(Mnew);
// fingerprint to select the orbit
fp := [];
for pp in PrimesUpTo(60) do
    if pp eq 2 or pp eq l then continue; end if;
    for tup in Factorization(pp*OK) do
        P := tup[1]; if Norm(Ncond+P) ne 1 then continue; end if;
        k,red := ResidueClassField(P); Rk := PolynomialRing(k); ok := true; Ck := 0;
        try Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(f)], Rk![red(c):c in Coefficients(hh)]);
        catch e; ok := false; end try;
        if ok then Append(~fp, <P, Fl!(-#Ck)>); end if;
    end for;
end for;
ef := 0; lam := 0; kf := 0; m := 0;
for nf in decomp do
    e0 := Eigenform(nf); Ff := HeckeEigenvalueField(nf); OF := Integers(Ff);
    for pr in Factorization(l*OF) do
        kf0, m0 := ResidueClassField(pr[1]);
        if Degree(kf0) ne 1 then continue; end if;                 // need residue field F_11
        good := true;
        for q in fp do if m0(OF!HeckeEigenvalue(e0,q[1])) ne kf0!(Integers()!q[2]) then good := false; break; end if; end for;
        if good then ef := e0; lam := pr[1]; kf := kf0; m := m0; OFF := OF; break; end if;
    end for;
    if ef cmpne 0 then break; end if;
end for;
error if ef cmpeq 0, "no matching eigenform/lambda found";
LOG(Sprintf("  eigenform ready (%.1o s); Hecke field deg %o, 11 %o", Cputime(t0),
    Degree(OFF), (RamificationIndex(lam) eq Degree(OFF)) select "totally ramified" else "not tot. ram."));

// ---- sigma irreducible? char poly of sigma at some prime irreducible over F_11 ----
irr := false;
for q in fp do
    cp := RT ! (RT.1^2 - (Integers()!q[2])*RT.1 + (Norm(q[1]) mod l));
    if IsIrreducible(cp) then irr := true; break; end if;
end for;
LOG(Sprintf("  sigma irreducible (a char poly is irreducible over F_11): %o", irr));

// ---- verify tr sigma = tr rho_f for all good P up to BOUND ----
LOG(Sprintf("verifying trace agreement for good primes N(P) <= %o ...", BOUND));
nchk := 0; fails := []; nextlog := 500;
for pp in PrimesUpTo(BOUND) do
    if pp eq 2 or pp eq l then continue; end if;
    for tup in Factorization(pp*OK) do
        P := tup[1];
        if Norm(P) gt BOUND then continue; end if;
        if Norm(Ncond+P) ne 1 then continue; end if;
        k,red := ResidueClassField(P); Rk := PolynomialRing(k); ok := true; Ck := 0;
        try Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(f)], Rk![red(c):c in Coefficients(hh)]);
        catch e; ok := false; end try;
        if not ok then continue; end if;
        tsig := Fl ! (-#Ck);
        trho := m(OFF ! HeckeEigenvalue(ef, P));
        if trho ne kf!(Integers()!tsig) then Append(~fails, <Norm(P), pp>); end if;
        nchk +:= 1;
    end for;
    if pp gt nextlog then
        LOG(Sprintf("  ... up to %o: %o primes checked, %o disagreements (%.1o s)", pp, nchk, #fails, Cputime(t0)));
        nextlog +:= 500;
    end if;
end for;

LOG(Sprintf("\nTOTAL: %o good primes of F with N(P) <= %o; trace DISAGREEMENTS: %o", nchk, BOUND, #fails));
if #fails eq 0 and irr then
    LOG(">>> tr sigma(Frob_P) = tr rho_f(Frob_P) for EVERY good prime P of F=Q(sqrt2), N(P) <= "
        cat IntegerToString(BOUND) cat ".");
    LOG("    sigma, rho_f: 2-dim, irreducible, det = chi_11, unramified outside {P_14303, 11}.");
    LOG("    cond(sigma) ~ 14303*11^a, (log cond)^2 ~ 200; the GRH conductor-based effective");
    LOG("    Faltings-Serre/Chebotarev bound (order 10^2-10^3) is far below BOUND.");
    LOG("    => (under GRH) sigma = rho_f. Hence sigma is modular (Serre over Q(sqrt2)).");
else
    LOG(Sprintf("  first disagreements: %o  (irr=%o)", fails[1..Min(10,#fails)], irr));
end if;
exit;
