// match_hmf_mod11.m -- Serre-modularity match, done ENTIRELY mod 11 (no number fields).
//
// At level = Curve A's bad prime above 66179, the cuspidal Hilbert modular forms space M
// (parallel weight 2) has integral Hecke operators T_P (Brandt matrices). A newform f with
// rho-bar_{f,lambda} = rho_i contributes, in M mod 11, an eigenvector v with T_P v = t_P v,
// where t_P = tr rho_i(Frob_P) is a root of  q_P(y) = y^2 - targetTr_P y + targetNm_P.
// Hence q_P(T_P) v = 0 for every test prime P, so
//     v  in  V := intersection over P of  ker( T_P^2 - targetTr_P T_P + targetNm_P I )  (mod 11).
// dim V > 0  <=>  rho_1/rho_2 ARE modular at this level  =>  evidence for Serre over Q(sqrt5).
// (targetTr,targetNm) from mod11_trnm.txt.
//
// Run: magma match_hmf_mod11.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;

K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
load "mod11_trnm.txt";   // trnm := [ <Np,p,rid,Tr,Nm>, ... ]

plog := "match_mod11_progress.txt";
PrintFile(plog, "# mod-11 HMF match progress" : Overwrite := true);
procedure LOG(s) PrintFile(plog, s); end procedure;

QI_A := [K|
    1, 1/1459240*(243125*w - 482787),
    1/2229718720*(-54798934*w + 127710391),
    1/8517525510400*(182422196780*w - 406668692089),
    1/65073894899456000*(38134182372761*w - 85270624049895) ];
bad := BadPrimesFromInvariants(QI_A, K);
Pbad := [P : P in bad | Norm(P) eq 66179][1];
p2 := Factorization(2*OK)[1][1];

function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;

// test primes: use a spread of small-norm primes first (cheap Hecke ops), coprime to level
targets := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    Append(~targets, <PrimeFromTag(row[2],row[3]), F11!row[4], F11!row[5], row[1], row[2]>);
end for;
Sort(~targets, func<a,b | a[4]-b[4]>);   // by norm ascending (cheapest Hecke ops first)

procedure CheckLevel(Level)
    printf "=== level norm %o ===\n", Norm(Level);
    LOG(Sprintf("level norm %o: building...", Norm(Level)));
    M := HilbertCuspForms(K, Level, [2,2]);
    d := Dimension(M);
    printf "  cuspidal dim = %o\n", d;
    LOG(Sprintf("  dim %o; intersecting kernels...", d));
    Vsp := VectorSpace(F11, d);    // current candidate space (start = everything)
    cur := Vsp;
    used := 0;
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;   // coprime to level
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        q := T*T - tg[2]*T + tg[3]*IdentityMatrix(F11, d);
        ker := Kernel(q);                      // subspace of Vsp
        cur := cur meet ker;
        used +:= 1;
        LOG(Sprintf("    after prime N=%o (p=%o): dim(intersection) = %o", tg[4], tg[5], Dimension(cur)));
        if Dimension(cur) eq 0 then
            printf "  intersection collapsed to 0 after %o primes -> NO match at this level.\n", used;
            LOG("    -> collapsed to 0, no match");
            return;
        end if;
        // once stable and small, enough primes to be conclusive
        if used ge 12 and Dimension(cur) le 4 then
            printf "  intersection STABLE at dim %o after %o primes.\n", Dimension(cur), used;
            break;
        end if;
    end for;
    printf "  *** dim V = %o after %o primes ***\n", Dimension(cur), used;
    if Dimension(cur) gt 0 then
        printf "  => rho_1/rho_2 ARE modular at level norm %o. SERRE EVIDENCE FOUND.\n", Norm(Level);
        LOG(Sprintf("  *** MATCH: dim V = %o at level norm %o ***", Dimension(cur), Norm(Level)));
    end if;
end procedure;

for Lev in [Pbad, Pbad*p2, Pbad*p2^2] do
    CheckLevel(Lev);
    printf "\n";
end for;
