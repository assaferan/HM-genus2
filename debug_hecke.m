// debug_hecke.m -- check the Hecke-eigenvalue / fingerprint convention at p=29.
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
load "mod11_trnm.txt";

QI_A := [K| 1, 1/1459240*(243125*w - 482787),
    1/2229718720*(-54798934*w + 127710391),
    1/8517525510400*(182422196780*w - 406668692089),
    1/65073894899456000*(38134182372761*w - 85270624049895) ];
Pbad := [P : P in BadPrimesFromInvariants(QI_A, K) | Norm(P) eq 66179][1];
p2 := Factorization(2*OK)[1][1];

function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;

// targets for p=29 (and a couple more small primes for cross-check)
function TargetsForP(pp)
    res := [];
    for row in trnm do
        if row[2] ne pp then continue; end if;
        P := PrimeFromTag(row[2], row[3]);
        Append(~res, <P, F11!row[4], F11!row[5], row[3]>);
    end for;
    return res;
end function;

for Lev in [Pbad, Pbad*p2] do
    printf "\n############ level norm %o ############\n", Norm(Lev);
    M := HilbertCuspForms(K, Lev, [2,2]);
    printf "dim = %o\n", Dimension(M);
    for pp in [29, 31, 41] do
        printf "\n--- p = %o ---\n", pp;
        for tg in TargetsForP(pp) do
            if Norm(tg[1] + Lev) ne 1 then printf "  (prime divides level, skip)\n"; continue; end if;
            T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
            evs := Eigenvalues(T);   // set of <eigenvalue, multiplicity> in F_11
            evset := { e[1] : e in evs };
            tTr := tg[2]; tNm := tg[3];
            R<y> := PolynomialRing(F11);
            qroots := { r[1] : r in Roots(y^2 - tTr*y + tNm) };
            printf "  prime rid=%o: targetTr=%o targetNm=%o  q-roots=%o\n",
                tg[4], tTr, tNm, qroots;
            printf "    q-roots in T-eigenvalues? %o ; -q-roots in? %o\n",
                qroots subset evset, {-r : r in qroots} subset evset;
            printf "    #distinct T-eigenvalues in F_11 = %o (of dim %o)\n", #evset, Dimension(M);
        end for;
    end for;
end for;
