// hmf_matching3.m -- Match rho_0 (2-dim mod-11 piece of Curve A) against
// NON-PARALLEL weight Hilbert modular forms over K = Q(sqrt(5)).
//
// Determinant constraint: det rho_0(Frob_p) = N(p) mod 11, while a weight
// (k1,k2) Hilbert form (trivial character) has det rho_f(Frob_p) = N(p)^{k0-1},
// k0 = max(k1,k2).  Matching for all good p forces k0 = 2 mod 10.  The smallest
// NON-parallel weight is therefore k0 = 12: weights [k1,12], k1 in {2,4,..,12}.
// (Check: N(p)^{11} = N(p) mod 11 by Fermat, so the determinant matches.)
//
// Filter (corrected): the Hecke eigenvalue is the TRACE a_p = alpha_1+alpha_2,
// so f lies in  cap_p ker(T_p - a_p I)  mod 11.  Trace values are weight-
// independent; we use the unambiguous INERT good primes.
//
// Run: magma hmf_matching3.m

K<s5> := QuadraticField(5);
OK    := Integers(K);
F11   := GF(11);

printf "=== rho_0 vs NON-PARALLEL weight HMF / Q(sqrt(5)) mod 11 ===\n\n";

// Inert good primes: <rational p, trace a_p mod 11>.
// a_p = -(middle coeff of char poly of Frob on rho_0) = alpha_1 + alpha_2.
//   p=7  : g=T^2+7T+5  -> a=4
//   p=13 : g=T^2+ T+4  -> a=10
//   p=37 : g=T^2   +5  -> a=0
//   p=47 : g=T^2+8T+9  -> a=3
//   p=53 : g=T^2+ T+4  -> a=10
inert := [ <7,4>, <13,10>, <37,0>, <47,3>, <53,10> ];

// Precompute inert prime ideals.
pid := [];
for e in inert do
    f := Factorization(e[1]*OK); assert #f eq 1;
    Append(~pid, <f[1][1], e[2]>);   // <prime ideal, a_p>
end for;

procedure TryWL(K, F11, pid, weight, level, label)
    printf "---------- weight %o, level %o (N=%o) ----------\n",
        weight, label, Norm(level);
    t0 := Cputime();
    M := HilbertCuspForms(K, level, weight);
    d := Dimension(M);
    printf "  dim = %o   (built in %.1os)\n", d, Cputime(t0);
    if d eq 0 then printf "  (empty)\n\n"; return; end if;

    V := VectorSpace(F11, d);
    cur := V;
    I := IdentityMatrix(F11, d);
    for pe in pid do
        P := pe[1]; ap := F11!pe[2];
        if Norm(level) mod Norm(P) eq 0 then continue; end if;
        Tp := ChangeRing(Matrix(HeckeOperator(M, P)), F11);
        // also report determinant of the 2-dim ad-hoc check: det piece = N(p)?
        ker := Kernel(hom<V->V | [V.i*(Tp - ap*I) : i in [1..d]]>);
        cur := cur meet ker;
        printf "    p=%o (a_p=%o): ker dim=%o, cumulative=%o\n",
            Norm(P), pe[2], Dimension(ker), Dimension(cur);
        if Dimension(cur) eq 0 then break; end if;
    end for;
    printf "  >>> surviving dim = %o\n\n", Dimension(cur);
end procedure;

// ----- candidate weights and levels -----
weights := [ [2,12], [12,2], [4,12], [12,4], [6,12], [12,6],
             [8,12], [12,8], [10,12], [12,10] ];

two := 2*OK;
fac := Factorization(66179*OK);
P1 := fac[1][1]; P2 := fac[2][1];

// Stage 1: small levels (fast) -- tests whether rho_0 is unramified at 66179.
printf "########## STAGE 1: small levels ##########\n\n";
for w in weights do
    TryWL(K, F11, pid, w, 1*OK, "(1)");
end for;
for w in weights do
    TryWL(K, F11, pid, w, two, "(2)");
end for;

// Stage 2: level 66179 (slow).
printf "########## STAGE 2: level 66179 ##########\n\n";
for w in weights do
    TryWL(K, F11, pid, w, P1, "P1[66179]");
end for;
for w in weights do
    TryWL(K, F11, pid, w, P2, "P2[66179]");
end for;

printf "=== done ===\n";
