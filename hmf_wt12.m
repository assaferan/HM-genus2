// hmf_wt12.m -- single (weight, level-prime) test for rho_0 vs non-parallel
// weight Hilbert forms over K = Q(sqrt(5)), level 66179.
//
// Invoke:  magma w1:=2 w2:=12 pidx:=1 hmf_wt12.m
//
// Filter: cap_p ker(T_p - a_p I) mod 11 over the residue field, INERT good
// primes, via INCREMENTAL subspace restriction (cheap after the first prime).
// Hecke matrices are over a fixed degree-4 number field; reduced entrywise
// modulo each prime above 11 (entries are integral in practice).

w1 := StringToInteger(w1); w2 := StringToInteger(w2);
pidx := StringToInteger(pidx);
weight := [w1, w2];

K<s5> := QuadraticField(5);
OK := Integers(K);
fac := Factorization(66179*OK);
level := fac[pidx][1];

// inert good primes with trace a_p = alpha_1 + alpha_2 (mod 11)
inert := [ <7,4>, <13,10>, <37,0>, <47,3>, <53,10> ];
pid := [];
for e in inert do
    fp := Factorization(e[1]*OK); assert #fp eq 1;
    Append(~pid, <fp[1][1], e[2]>);
end for;

printf "=== weight %o, level P%o[66179] ===\n", weight, pidx;
t0 := Cputime();
M := HilbertCuspForms(K, level, weight);
d := Dimension(M);
printf "dim = %o  (built %.1os)\n", d, Cputime(t0); ttt := Cputime();
if d eq 0 then printf "empty\n"; exit; end if;

// integral entrywise reduction of matrix over R modulo Q (residue map mp:OR->k)
function ReduceAt(TM, k, mp, OR)
    n := Nrows(TM);
    return Matrix(k, n, n, [mp(OR!x) : x in Eltseq(TM)]);
end function;

// compute & cache raw Hecke matrices over R once (reused for both Q|11)
printf "Computing %o Hecke operators over the base ring ...\n", #pid;
raw := [* *];
for pe in pid do
    th := Cputime();
    Append(~raw, Matrix(HeckeOperator(M, pe[1])));
    printf "  HeckeOp p=%o : %.1os\n", Norm(pe[1]), Cputime(th);
end for;
R := BaseRing(raw[1]);
printf "Hecke base ring: %o\n", R;
OR := Integers(R);
fac11 := Factorization(11*OR);
printf "11 -> %o prime(s), residue degrees %o\n",
    #fac11, [Degree(ResidueClassField(q[1])) : q in fac11];

// For each prime Q|11: incremental restriction.
for qf in fac11 do
    Q := qf[1];
    k, mp := ResidueClassField(Q);
    printf "-- Q|11, residue field %o --\n", k;
    n := d;
    W := IdentityMatrix(k, n);  // current surviving basis (rows), full to start
    Wset := false;
    for i in [1..#pid] do
        pe := pid[i]; ap := k ! pe[2];
        tr := Cputime();
        A := ReduceAt(raw[i], k, mp, OR);
        Am := A - ScalarMatrix(k, n, ap);
        if not Wset then
            ker := Kernel(Am);              // first prime: full n x n kernel
            W := BasisMatrix(ker);
            Wset := true;
        else
            C := W * Am;                    // (w x n)
            ker := Kernel(C);               // x in k^w with x*C = 0
            W := BasisMatrix(ker) * W;
        end if;
        printf "  p=%o (a_p=%o): surviving dim=%o   (%.1os)\n",
            Norm(pe[1]), pe[2], Nrows(W), Cputime(tr);
        if Nrows(W) eq 0 then break; end if;
    end for;
    printf "  >>> MATCHED dim over %o = %o\n", k, Nrows(W);
end for;
printf "=== done weight %o P%o (total %.1os) ===\n", weight, pidx, Cputime(ttt);
