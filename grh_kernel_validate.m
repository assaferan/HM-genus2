// grh_kernel_validate.m -- Step 1: validate the "survivor eigenvector" route for computing
// a_P(f) mod lambda WITHOUT isolating the char-0 eigenform.
//
// Claim: the 1-dim kernel survivor v (from kernel_torsion.m) is a simultaneous mod-l Hecke
// eigenvector, and for every good prime Q the scalar c_Q with  v * T_Q = c_Q * v  equals
// a_Q(f) mod lambda -- i.e. exactly the quantity the GRH certificate (grh_14303.m) currently
// gets from the char-0 eigenform.  Here we PROVE it on 14303.1, where we DO have the char-0
// eigenform, by checking  c_Q  ==  reduce(HeckeEigenvalue(ef,Q))  ==  -#C(F_Q)  for all good
// Q up to BOUND.  If they agree, the survivor route is trustworthy where no char-0 form exists.
//
// Success sentinel:  GRH-KERNEL VALIDATE: PASS
// Run: magma grh_kernel_validate.m           (optional BOUND:=1000 to go faster)
SetColumns(0);
if not assigned BOUND then BOUND := 3000; end if;
BOUND := StringToInteger(Sprintf("%o", BOUND));
K<a> := QuadraticField(2); OK := Integers(K); R<x> := PolynomialRing(K);
l := 11; Fl := GF(l);

// curve 14303.1
f := a*x^5 + (3+2*a)*x^4 + (3+a)*x^3 + (-4-2*a)*x^2 + (-1+2*a)*x + (1-a);
hh := x^3 + x^2 + 1;
C := HyperellipticCurve(f, hh);
N := ideal<OK | -13*a - 121>;                     // norm 14303, e=0
assert Norm(N) eq 14303;

// fingerprint t_P = -#C(F_P) mod l
function tP(P)
    k,red := ResidueClassField(P); Rk := PolynomialRing(k);
    Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(f)], Rk![red(c):c in Coefficients(hh)]);
    return Fl!(-#Ck);
end function;
fpP := []; fpT := [];
for pp in PrimesUpTo(60) do
    if #fpP ge 12 then break; end if;
    if pp eq 2 or pp eq l then continue; end if;
    for tup in Factorization(pp*OK) do
        P := tup[1]; if Norm(N+P) ne 1 then continue; end if;
        Append(~fpP, P); Append(~fpT, tP(P));
        if #fpP ge 12 then break; end if;
    end for;
end for;

// ---- survivor eigenvector v on the FULL space (the kernel-method object) ----
printf "building full space + survivor ...\n"; t0 := Cputime();
M := HilbertCuspForms(K, N, [2,2]); dm := Dimension(M);
assert dm eq 595;
V := VectorSpace(Fl, dm); cur := V; I := IdentityMatrix(Fl, dm);
for i in [1..#fpP] do
    Tp := ChangeRing(Matrix(HeckeOperator(M, fpP[i])), Fl);
    cur := cur meet Kernel(Tp - fpT[i]*I);
end for;
assert Dimension(cur) eq 1;                       // isolates a single eigensystem
v := Basis(cur)[1];
j := rep{ i : i in [1..dm] | v[i] ne 0 };         // pivot coordinate for reading eigenvalues
printf "  survivor is 1-dim (%.1o s)\n", Cputime(t0);

// eigenvalue of T_Q on the survivor:  v * T_Q = c * v
function survEig(Q)
    Tq := ChangeRing(Matrix(HeckeOperator(M, Q)), Fl);
    w := v*Tq;
    c := w[j]/v[j];
    error if w ne c*v, Sprintf("survivor NOT an eigenvector at N(Q)=%o", Norm(Q));
    return c;
end function;

// ---- char-0 eigenform (ground truth), selected by the fingerprint ----
printf "building char-0 eigenform (ground truth) ...\n";
Mnew := NewSubspace(M); decomp := NewformDecomposition(Mnew);
ef := 0; lam := 0; kf := 0; m := 0; OFF := 0;
for nf in decomp do
    e0 := Eigenform(nf); Ff := HeckeEigenvalueField(nf); OF := Integers(Ff);
    for pr in Factorization(l*OF) do
        kf0, m0 := ResidueClassField(pr[1]);
        if Degree(kf0) ne 1 then continue; end if;
        good := true;
        for i in [1..#fpP] do
            if m0(OF!HeckeEigenvalue(e0,fpP[i])) ne kf0!(Integers()!fpT[i]) then good := false; break; end if;
        end for;
        if good then ef := e0; lam := pr[1]; kf := kf0; m := m0; OFF := OF; break; end if;
    end for;
    if ef cmpne 0 then break; end if;
end for;
error if ef cmpeq 0, "no matching char-0 eigenform";
printf "  char-0 eigenform ready; Hecke field deg %o\n", Degree(OFF);

// ---- compare survivor-eigenvalue vs char-0 reduction vs curve fingerprint ----
printf "comparing all three at good primes N(P) <= %o ...\n", BOUND;
nchk := 0; bad := [];
for pp in PrimesUpTo(BOUND) do
    if pp eq 2 or pp eq l then continue; end if;
    for tup in Factorization(pp*OK) do
        P := tup[1];
        if Norm(P) gt BOUND or Norm(N+P) ne 1 then continue; end if;
        cv   := survEig(P);                        // survivor route (no char-0 form)
        c0   := m(OFF ! HeckeEigenvalue(ef, P));   // char-0 reduction (ground truth)
        tsig := kf!(Integers()!tP(P));             // curve fingerprint
        if not (cv eq c0 and c0 eq tsig) then Append(~bad, <Norm(P), cv, c0, tsig>); end if;
        nchk +:= 1;
    end for;
end for;

printf "\nchecked %o good primes; mismatches: %o\n", nchk, #bad;
if #bad eq 0 then
    printf "GRH-KERNEL VALIDATE: PASS (survivor-eigenvalue = char-0 reduction = fingerprint at all %o primes, N(P)<=%o)\n", nchk, BOUND;
else
    printf "  first mismatches <N(P),surv,char0,fp>: %o\n", bad[1..Min(10,#bad)];
    printf "GRH-KERNEL VALIDATE: FAIL\n";
end if;
exit;
