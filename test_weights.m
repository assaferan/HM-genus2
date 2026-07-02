// test_weights.m -- search for the genuine modular form realizing rho_i at NON-PARALLEL
// weights (11 splits => Serre weight can differ at the two primes above 11).
//
// For each weight, V+ = cap_P ker(T_P^2 - tTr*T_P + tNm) (untwisted: a_P mod 11 in {t1,t2}),
// Vtw = cap_P ker[(T^2 - tTr T + tNm)(T^2 + tTr T + tNm)] (allow quadratic twist).
// dim V+ > 0  =>  candidate genuine match at this weight.
//
// Run: magma test_weights.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
load "mod11_trnm.txt";

LOGF := "test_weights_out.txt";
PrintFile(LOGF, "# test_weights: non-parallel weight search" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI_A := [K| 1, 1/1459240*(243125*w - 482787),
    1/2229718720*(-54798934*w + 127710391),
    1/8517525510400*(182422196780*w - 406668692089),
    1/65073894899456000*(38134182372761*w - 85270624049895) ];
Pbad := [P : P in BadPrimesFromInvariants(QI_A, K) | Norm(P) eq 66179][1];
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

procedure TryWeight(wt, NFIND)
    LOG(Sprintf("==== weight %o, level Pbad ====", wt));
    t0 := Cputime();
    M := HilbertCuspForms(K, Pbad, wt);
    d := Dimension(M);
    LOG(Sprintf("  dim = %o (build %o s)", d, Cputime(t0)));
    if d eq 0 then LOG("  empty."); return; end if;
    if d gt 40000 then LOG("  too big, skipping Hecke."); return; end if;
    Id := IdentityMatrix(F11, d);
    Vp := VectorSpace(F11, d); Vt := Vp; nf := 0;
    for tg in targets do
        if Norm(tg[1] + Pbad) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11); T2 := T*T;
        kp := Kernel(T2 - tg[2]*T + tg[3]*Id);
        km := Kernel(T2 + tg[2]*T + tg[3]*Id);
        Vp := Vp meet kp; Vt := Vt meet (kp + km);
        T := 0; T2 := 0; kp := 0; km := 0; nf +:= 1;
        LOG(Sprintf("    N=%o: dim V+=%o  Vtw=%o", tg[4], Dimension(Vp), Dimension(Vt)));
        if Dimension(Vt) eq 0 then break; end if;
        if nf ge NFIND then break; end if;
    end for;
    LOG(Sprintf("  => weight %o: dim V+ = %o, dim Vtw = %o", wt, Dimension(Vp), Dimension(Vt)));
    if Dimension(Vp) gt 0 then LOG("  *** UNTWISTED MATCH at this weight! ***"); end if;
end procedure;

TryWeight([2,12], 12);
TryWeight([12,2], 12);
LOG("DONE.");
exit;
