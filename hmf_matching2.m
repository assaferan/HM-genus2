// hmf_matching2.m -- Match the 2-dim mod-11 Galois rep rho_0 of Curve A
// against Hilbert modular forms over K = Q(sqrt(5)), the RIGHT way:
//
//   f lies in  intersection_p  ker( g_p(T_p) )   (mod 11),
//
// where g_p is the irreducible quadratic = char poly of Frob_p on rho_0.
// No traces, no eigenvalue extraction: we plug the matrix T_p straight into
// the quadratic and take the kernel.
//
// We restrict to INERT good primes {7,13,37,47,53} so the prime-ideal <-> g_p
// matching is unambiguous (one prime above p, norm p^2).
//
// Since cond(rho_0) away from 11 is unknown (it divides 2^a * 66179^b, the odd
// part of cond(A) being supported only at 66179), we test several levels.
//
// Run: magma hmf_matching2.m

K<s5> := QuadraticField(5);
OK    := Integers(K);
F11   := GF(11);
PZ<T> := PolynomialRing(F11);

printf "=== HMF matching for rho_0 (Curve A) / Q(sqrt(5)) mod 11 ===\n";
printf "    using f in  cap_p ker(g_p(T_p)),  g_p = irreducible quadratic.\n\n";

// -----------------------------------------------------------------------
// L-polynomials L_p(T) of Curve A at the inert good primes, as coefficient
// lists [c0, c1, c2, c3, c4] with L_p(T) = sum c_k T^k  (c0 = 1).
// (Transcribed from data_qsqrt5_gp11.txt; recomputed factors below.)
// -----------------------------------------------------------------------
inert_data := [
    <7,  [1,   2,  -58,    98,    2401   ]>,
    <13, [1,  17,  376,    2873,  28561  ]>,
    <37, [1,  -6,  -1310, -8214,  1874161]>,
    <47, [1,   9,  1808,   19881, 4879681]>,
    <53, [1,  50,  2950,   140450,7890481]>
];

// For each inert prime build:
//   g_frob = irreducible quadratic factor of the CHAR POLY OF FROBENIUS
//            (reciprocal of L_p), roots = Frobenius eigenvalues alpha_i.
//   g_lpol = irreducible quadratic factor of L_p itself, roots = 1/alpha_i.
// Both as monic quadratics over F11.
prime_info := [];  // list of <p, prime_ideal, g_frob, g_lpol>
printf "Quadratic factors g_p (both normalizations):\n";
for entry in inert_data do
    p := entry[1];
    c := entry[2];
    // char poly of Frobenius = X^4 + c1 X^3 + c2 X^2 + c3 X + c4 (monic)
    cp_frob := T^4 + (F11!c[2])*T^3 + (F11!c[3])*T^2 + (F11!c[4])*T + (F11!c[5]);
    // L_p mod 11 (monic-ized)
    Lp := (F11!c[5])*T^4 + (F11!c[4])*T^3 + (F11!c[3])*T^2 + (F11!c[2])*T + (F11!c[1]);
    Lp := Lp / LeadingCoefficient(Lp);

    g_frob := [f[1] : f in Factorization(cp_frob) | Degree(f[1]) eq 2][1];
    g_lpol := [f[1] : f in Factorization(Lp)      | Degree(f[1]) eq 2][1];

    fac_p := Factorization(p*OK);
    assert #fac_p eq 1;          // inert
    pid := fac_p[1][1];

    Append(~prime_info, <p, pid, g_frob, g_lpol>);
    printf "  p=%o (N=%o):  g_frob = %o ,  g_lpol = %o\n",
        p, Norm(pid), g_frob, g_lpol;
end for;
printf "\n";

// -----------------------------------------------------------------------
// For a candidate level ideal: build S_2(K, level), then intersect
//   ker(g_p(T_p))  over the inert good primes, separately for the two
//   normalizations of g_p.  Skip primes dividing the level.
// -----------------------------------------------------------------------
procedure TryLevel(K, F11, prime_info, level, label)
    printf "================ Level %o  (N = %o) ================\n", label, Norm(level);
    M := HilbertCuspForms(K, level);
    d := Dimension(M);
    printf "dim S_2(K, level) = %o\n", d;
    if d eq 0 then printf "(empty space)\n\n"; return; end if;

    V := VectorSpace(F11, d);
    cur_frob := V;
    cur_lpol := V;
    I := IdentityMatrix(F11, d);

    for info in prime_info do
        p := info[1]; pid := info[2]; gf := info[3]; gl := info[4];
        if Norm(level) mod p eq 0 then
            printf "  [skip p=%o : divides level]\n", p; continue;
        end if;
        Tp := ChangeRing(Matrix(HeckeOperator(M, pid)), F11);

        // g(T_p) = Tp^2 + u Tp + w I  for g = T^2 + u T + w
        uf := Coefficient(gf,1); wf := Coefficient(gf,0);
        ul := Coefficient(gl,1); wl := Coefficient(gl,0);
        Mf := Tp*Tp + uf*Tp + wf*I;
        Ml := Tp*Tp + ul*Tp + wl*I;

        kf := Kernel(hom<V->V | [V.i*Mf : i in [1..d]]>);
        kl := Kernel(hom<V->V | [V.i*Ml : i in [1..d]]>);
        cur_frob := cur_frob meet kf;
        cur_lpol := cur_lpol meet kl;

        printf "  p=%o:  dim ker g_frob(T_p)=%o (cum %o) | dim ker g_lpol(T_p)=%o (cum %o)\n",
            p, Dimension(kf), Dimension(cur_frob), Dimension(kl), Dimension(cur_lpol);
    end for;

    printf ">>> Level %o : final surviving dim  [g_frob]=%o   [g_lpol]=%o\n\n",
        label, Dimension(cur_frob), Dimension(cur_lpol);
end procedure;

// -----------------------------------------------------------------------
// Candidate levels.  2 is INERT in Q(sqrt(5)) (norm 4).  66179 SPLITS.
// -----------------------------------------------------------------------
two := 2*OK;                                   // inert prime above 2, N=4
fac66179 := Factorization(66179*OK);
P1 := fac66179[1][1];
P2 := fac66179[2][1];

// Small levels first (fast), then the 66179 ones (slow ~minutes each).
TryLevel(K, F11, prime_info, 1*OK,        "(1)");
TryLevel(K, F11, prime_info, two,         "(2)");
TryLevel(K, F11, prime_info, two^2,       "(2)^2");
TryLevel(K, F11, prime_info, two^3,       "(2)^3");
TryLevel(K, F11, prime_info, P1,          "P1[66179]");
TryLevel(K, F11, prime_info, P2,          "P2[66179]");
TryLevel(K, F11, prime_info, two*P1,      "(2)*P1");
TryLevel(K, F11, prime_info, two*P2,      "(2)*P2");
TryLevel(K, F11, prime_info, P1*P2,       "P1*P2=(66179)");

printf "=== done ===\n";
