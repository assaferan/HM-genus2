// hmf_matching.m -- Match the 2-dimensional mod-11 Galois representation of
// Curve A (from the (1,11) Gross-Popescu search over Q(sqrt(5))) against
// Hilbert modular forms.
//
// Background: by construction the mod-11 representation on Jac[11] splits
// as 2 + 1 + 1.  At any prime P where the 4-dim L-polynomial mod 11 has an
// irreducible quadratic factor, that quadratic is the L-polynomial of the 2-dim
// piece (the other two factors are linear).
//
// Key conversion: if the irreducible quadratic factor in the *monic* L_P(T) mod 11
// is t^2 + a*t + b, then the Frobenius eigenvalues alpha_1, alpha_2 on the 2-dim
// piece satisfy:
//   1/alpha_i are the roots of t^2+at+b
//   => alpha_1+alpha_2 = -a/b  (mod 11)    [= Hecke eigenvalue a_P]
//   => alpha_1*alpha_2 = 1/b = N(P) mod 11  [check: matches cyclotomic det]
//
// So the filter on the HMF space is the EIGENSPACE of T_P for eigenvalue a_P = -a/b.
//
// Run: magma hmf_matching.m

AttachSpec("../CHIMP/CHIMP.spec");
load "survivors_qsqrt5_gp11.m";
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, ConductorExponentAt;

K<s5> := QuadraticField(5);
OK := Integers(K);
F11 := GF(11);

printf "=== HMF matching for Curve A / Q(sqrt(5)) mod 11 ===\n\n";

// -----------------------------------------------------------------------
// Step 1.  Identify both prime ideals above 66179.
// 66179 splits in Q(sqrt(5)) because (5/66179) = (66179 mod 5 / 5) = (4/5) = 1.
// We try BOTH primes as the HMF level.
// -----------------------------------------------------------------------
fac66179 := Factorization(66179 * OK);
assert #fac66179 eq 2 and &and[f[2] eq 1 : f in fac66179];
P1 := fac66179[1][1];
P2 := fac66179[2][1];
printf "Primes above 66179: P1 (N=%o), P2 (N=%o)\n\n", Norm(P1), Norm(P2);

// -----------------------------------------------------------------------
// Steps 2-5 run for each candidate level prime.
// -----------------------------------------------------------------------

procedure TryLevel(K, OK, F11, cond_N, constraints)
    printf "=== Trying level N = prime above 66179, N(N) = %o ===\n", Norm(cond_N);

    // Step 2.  Construct S_2(K, N) and report dimension.
    printf "Building S_2(K, N) ...\n";
    M := HilbertCuspForms(K, cond_N);
    d := Dimension(M);
    printf "dim S_2(K, N) = %o\n\n", d;
    if d eq 0 then printf "Space is zero -- no HMF at this level.\n\n"; return; end if;

    // Step 4.  Intersect eigenspaces of T_P mod 11 iteratively.
    V   := VectorSpace(F11, d);
    cur := V;
    printf "Intersecting eigenspace constraints mod 11 (initial dim = %o):\n", d;

    for entry in constraints do
        p_rat := entry[1]; idx := entry[2]; ap := F11 ! entry[3];
        fac_p := Factorization(p_rat * OK);
        if #fac_p lt idx then
            printf "  [skip p=%o idx=%o: only %o primes above p]\n", p_rat, idx, #fac_p;
            continue;
        end if;
        pid := fac_p[idx][1];

        T_Z  := Matrix(HeckeOperator(M, pid));
        T11  := ChangeRing(T_Z, F11);
        I    := IdentityMatrix(F11, d);

        ker := Kernel(hom<V -> V | [V.i * (T11 - ap*I) : i in [1..d]]>);
        cur := cur meet ker;

        printf "  T_{p=%o, N=%o, #%o}: a_P = %o  ->  surviving dim = %o\n",
            p_rat, Norm(pid), idx, entry[3], Dimension(cur);

        if Dimension(cur) eq 0 then
            printf "  !! Dimension 0 after this constraint.\n";
            // Print actual eigenvalues of the form selected only by the first (p=7) constraint.
            printf "  [Printing actual eigenvalues of form(s) with correct a_{P_7}=4:]\n";
            V2 := VectorSpace(F11, d);
            cur2 := V2;
            entry0 := constraints[1];
            fac0 := Factorization(entry0[1] * OK);
            pid0 := fac0[entry0[2]][1];
            T0 := ChangeRing(Matrix(HeckeOperator(M, pid0)), F11);
            ker0 := Kernel(hom<V2 -> V2 | [V2.i * (T0 - (F11!entry0[3])*IdentityMatrix(F11,d)) : i in [1..d]]>);
            if Dimension(ker0) gt 0 then
                v := Basis(ker0)[1];
                for e2 in constraints do
                    p2 := e2[1]; idx2 := e2[2];
                    fac2 := Factorization(p2 * OK);
                    if #fac2 lt idx2 then continue; end if;
                    pid2 := fac2[idx2][1];
                    T2 := ChangeRing(Matrix(HeckeOperator(M, pid2)), F11);
                    vT := v * T2;
                    lambda := F11!0;
                    for ci in [1..d] do
                        if v[ci] ne 0 then lambda := vT[ci] / v[ci]; break; end if;
                    end for;
                    printf "    T_{p=%o,N=%o,#%o}: actual=%o  predicted=%o  match:%o\n",
                        p2, Norm(pid2), idx2, Integers()!lambda, e2[3], lambda eq F11!e2[3];
                end for;
            end if;
            break;
        end if;

        if Dimension(cur) le 5 then
            printf "    [dim small -- actual eigenvalues:]\n";
            v := Basis(cur)[1];
            for e2 in constraints do
                p2 := e2[1]; idx2 := e2[2];
                fac2 := Factorization(p2 * OK);
                if #fac2 lt idx2 then continue; end if;
                pid2 := fac2[idx2][1];
                T2 := ChangeRing(Matrix(HeckeOperator(M, pid2)), F11);
                vT := v * T2;
                lambda := F11!0;
                for ci in [1..d] do
                    if v[ci] ne 0 then lambda := vT[ci] / v[ci]; break; end if;
                end for;
                printf "    T_{p=%o,N=%o,#%o}: actual=%o  predicted=%o  match:%o\n",
                    p2, Norm(pid2), idx2, Integers()!lambda, e2[3], lambda eq F11!e2[3];
            end for;
        end if;
    end for;

    printf "\n=== Final surviving dim (level N(N)=%o): %o ===\n\n",
        Norm(cond_N), Dimension(cur);

    if Dimension(cur) gt 0 then
        printf "Decomposing into newform spaces ...\n";
        NF := NewformDecomposition(M);
        printf "Total newform spaces: %o\n\n", #NF;
        good := [];
        for i in [1..#NF] do
            S := NF[i];
            ok := true;
            for entry in constraints do
                p_rat := entry[1]; idx := entry[2]; ap := entry[3];
                fac_p := Factorization(p_rat * OK);
                if #fac_p lt idx then continue; end if;
                pid := fac_p[idx][1];
                T_S := ChangeRing(Matrix(HeckeOperator(S, pid)), F11);
                dim_S := Nrows(T_S);
                if T_S ne (F11!ap) * IdentityMatrix(F11, dim_S) then ok := false; break; end if;
            end for;
            if ok then
                Append(~good, i);
                printf "  Newform space %o (dim %o): all constraints satisfied\n", i, Dimension(S);
            end if;
        end for;
        printf "\n%o newform space(s) match all constraints at this level.\n\n", #good;
    end if;
end procedure;

// -----------------------------------------------------------------------
// Step 3.  Hecke eigenvalues from Curve A Frobenius data.
//
// For each usable prime P (irreducible quadratic factor in L_P mod 11),
// the quadratic factor is t^2 + a*t + b in the MONIC L_P(T) mod 11.
// Hecke eigenvalue: a_P = -a * b^(-1) mod 11.
//
// Format: <rational prime p, prime index i, a_P mod 11>
// (We also record [a,b] for cross-check.)
// -----------------------------------------------------------------------
// Pre-computed from the Frobenius data in data_qsqrt5_gp11.txt:
//   N(P)=49  (p=7  inert): t^2+8t+9  -> -8*9^(-1) = -8*5 = -40 =  4 mod 11
//   N(P)=169 (p=13 inert): t^2+3t+3  -> -3*3^(-1) = -1         = 10 mod 11
//   N(p1)=29 (split, #1):  t^2+10t+8 -> -10*8^(-1)= -10*7=-70  =  7 mod 11
//   N(p2)=29 (split, #2):  t^2+7t+8  -> -7*8^(-1) = -7*7=-49   =  6 mod 11
//   N(p1)=31 (split, #1):  t^2+9t+5  -> -9*5^(-1) = -9*9=-81   =  7 mod 11
//   N(P)=1369(p=37 inert): t^2+0t+9  ->  0*9^(-1) =  0         =  0 mod 11
//   N(p2)=41 (split, #2):  t^2+5t+7  -> -5*7^(-1) = -5*8=-40   =  4 mod 11
//   N(P)=2209(p=47 inert): t^2+7t+5  -> -7*5^(-1) = -7*9=-63   =  3 mod 11
//   N(P)=2809(p=53 inert): t^2+3t+3  -> -3*3^(-1) = -1         = 10 mod 11
// Determinant check: alpha_1*alpha_2 = b^(-1) should equal N(P) mod 11 -- verified.
constraints := [
    <7,  1,  4>,   // N(P)=49,   a_P = 4   (7 inert)
    <13, 1, 10>,   // N(P)=169,  a_P = 10  (13 inert)
    <29, 1,  7>,   // N(p1)=29,  a_p1 = 7  (29 split, prime 1)
    <29, 2,  6>,   // N(p2)=29,  a_p2 = 6  (29 split, prime 2)
    <31, 1,  7>,   // N(p1)=31,  a_p1 = 7  (31 split, prime 1)
    <37, 1,  0>,   // N(P)=1369, a_P = 0   (37 inert)
    <41, 2,  4>,   // N(p2)=41,  a_p2 = 4  (41 split, prime 2)
    <47, 1,  3>,   // N(P)=2209, a_P = 3   (47 inert)
    <53, 1, 10>    // N(P)=2809, a_P = 10  (53 inert)
];

TryLevel(K, OK, F11, P1, constraints);
TryLevel(K, OK, F11, P2, constraints);
