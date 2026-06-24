// match_hmf_twist.m -- twist-robust Serre-modularity match, mod 11.
//
// The debug showed the matching form realizes rho_i UP TO A QUADRATIC TWIST eps (ramified
// at 2): a_P(f) = eps(P) * t_P with t_P a root of q_P^+ = y^2 - tTr y + tNm. So a_P(f) is a
// root of EITHER q_P^+ or q_P^- = y^2 + tTr y + tNm. Hence the eigenform lies in
//     ker q_P^+(T_P)  +  ker q_P^-(T_P)
// at every prime P, and the matching space is V = intersection_P of that sum.
// dim V > 0  =>  rho_i (up to quadratic twist) is modular  =>  Serre evidence.
//
// Run: magma match_hmf_twist.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
load "mod11_trnm.txt";

plog := "match_twist_progress.txt";
PrintFile(plog, "# twist-robust mod-11 match" : Overwrite := true);
procedure LOG(s) PrintFile(plog, s); end procedure;

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

targets := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    Append(~targets, <PrimeFromTag(row[2],row[3]), F11!row[4], F11!row[5], row[1], row[2]>);
end for;
Sort(~targets, func<a,b | a[4]-b[4]>);

procedure CheckLevel(Level, maxprimes)
    printf "=== level norm %o ===\n", Norm(Level);
    LOG(Sprintf("level norm %o: building...", Norm(Level)));
    M := HilbertCuspForms(K, Level, [2,2]);
    d := Dimension(M);
    printf "  cuspidal dim = %o\n", d;
    Id := IdentityMatrix(F11, d);
    cur := VectorSpace(F11, d); used := 0; stable := 0; lastdim := d;
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        T2 := T*T;
        kplus  := Kernel(T2 - tg[2]*T + tg[3]*Id);    // untwisted
        kminus := Kernel(T2 + tg[2]*T + tg[3]*Id);    // twisted by -1
        dkp := Dimension(kplus); dkm := Dimension(kminus);
        cur := cur meet (kplus + kminus);
        T := 0; T2 := 0; kplus := 0; kminus := 0;     // free ~30MB matrices each iter
        used +:= 1;
        if Dimension(cur) eq lastdim then stable +:= 1; else stable := 0; end if;
        lastdim := Dimension(cur);
        LOG(Sprintf("    N=%o (p=%o): dim+=%o dim-=%o  dim(V)=%o",
            tg[4], tg[5], dkp, dkm, Dimension(cur)));
        printf "  N=%o: dim(V) = %o\n", tg[4], Dimension(cur);
        if Dimension(cur) eq 0 then printf "  -> collapsed to 0.\n"; LOG("  collapsed"); return; end if;
        if used ge maxprimes or (stable ge 6 and Dimension(cur) le 4) then break; end if;
    end for;
    printf "  *** dim V = %o after %o primes ***\n", Dimension(cur), used;
    LOG(Sprintf("  *** dim V = %o (%o primes) ***", Dimension(cur), used));
    if Dimension(cur) gt 0 then
        printf "  => SERRE EVIDENCE: rho_i (up to quadratic twist at 2) IS modular at level norm %o.\n", Norm(Level);
    end if;
end procedure;

CheckLevel(Pbad*p2, 14);
