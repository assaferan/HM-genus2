// test_level2.m -- weight (2,2) is confirmed (ordinary at both primes|11). Search the genuine
// UNTWISTED match at higher 2-part levels Pbad*(2)^2 (and Pbad*(2)^3 if feasible).
//   V+ = cap_P ker(T_P^2 - tTr T_P + tNm)  (a_P mod 11 a root of q+ = tr rho_i)
// dim V+ > 0 => rho_i modular by a weight-2 form at this level.
//
// Run: magma test_level2.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
load "mod11_trnm.txt";

LOGF := "test_level2_out.txt";
PrintFile(LOGF, "# test_level2: untwisted match at higher 2-part level, weight (2,2)" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

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

procedure TryLevel(Level, NFIND)
    LOG(Sprintf("==== level norm %o ====", Norm(Level)));
    t0 := Cputime();
    M := HilbertCuspForms(K, Level, [2,2]);
    d := Dimension(M);
    LOG(Sprintf("  dim = %o (build %o s)", d, Cputime(t0)));
    if d eq 0 then return; end if;
    Id := IdentityMatrix(F11, d);
    Vp := VectorSpace(F11, d); Vt := Vp; nf := 0;
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        th := Cputime();
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11); T2 := T*T;
        kp := Kernel(T2 - tg[2]*T + tg[3]*Id);
        km := Kernel(T2 + tg[2]*T + tg[3]*Id);
        Vp := Vp meet kp; Vt := Vt meet (kp + km);
        T := 0; T2 := 0; kp := 0; km := 0; nf +:= 1;
        LOG(Sprintf("    N=%o: dim V+=%o Vtw=%o  (Hecke %o s)", tg[4], Dimension(Vp), Dimension(Vt), Cputime(th)));
        if Dimension(Vt) eq 0 then break; end if;
        if nf ge NFIND then break; end if;
    end for;
    LOG(Sprintf("  => level norm %o: dim V+ = %o, Vtw = %o", Norm(Level), Dimension(Vp), Dimension(Vt)));
    if Dimension(Vp) gt 0 then LOG("  *** UNTWISTED MATCH: rho_i modular, weight (2,2), this level! ***"); end if;
end procedure;

TryLevel(Pbad*p2^2, 12);
LOG("DONE.");
exit;
