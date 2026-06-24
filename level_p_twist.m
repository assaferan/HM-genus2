// level_p_twist.m -- twist-robust search at level p (norm 66179, dim 1102), looking for a
// form realizing rho_2 (or rho_1 untwisted). Cheap. Also extract a_P and label vs the
// fingerprint roots.
//
// Run: magma level_p_twist.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K); F11 := GF(11); R<y> := PolynomialRing(F11);
load "mod11_trnm.txt";

QI_A := [K| 1, 1/1459240*(243125*w - 482787),
    1/2229718720*(-54798934*w + 127710391),
    1/8517525510400*(182422196780*w - 406668692089),
    1/65073894899456000*(38134182372761*w - 85270624049895) ];
Pbad := [P : P in BadPrimesFromInvariants(QI_A, K) | Norm(P) eq 66179][1];
function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if; return ideal<OK | p, w - rid>; end function;
targets := [];
for r in trnm do
    if r[2] eq 11 then continue; end if;
    Append(~targets, <PrimeFromTag(r[2],r[3]), F11!r[4], F11!r[5], r[1], r[2], r[3]>);
end for;
Sort(~targets, func<a,b | a[4]-b[4]>);

// try BOTH primes above 66179 (Pbad = bad one; Pbad2 = the other, good for A)
Pall := [fa[1] : fa in Factorization(66179*OK)];
Level := [P : P in Pall | P ne Pbad][1];   // the OTHER prime above 66179
printf "level = OTHER prime above 66179, norm %o\n", Norm(Level);
M := HilbertCuspForms(K, Level, [2,2]); d := Dimension(M); Id := IdentityMatrix(F11, d);
printf "cuspidal dim = %o\n", d;
cur := VectorSpace(F11, d); built := 0;
for tg in targets do
    if Norm(tg[1] + Level) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11); T2 := T*T;
    kp := Kernel(T2 - tg[2]*T + tg[3]*Id); km := Kernel(T2 + tg[2]*T + tg[3]*Id);
    cur := cur meet (kp + km); T:=0; T2:=0; kp:=0; km:=0;
    built +:= 1;
    printf "  after N=%o: dim V = %o\n", tg[4], Dimension(cur);
    if Dimension(cur) eq 0 then printf "  -> collapsed; rho_2 NOT at level p.\n"; exit; end if;
    if built ge 10 then break; end if;
end for;
printf "\n*** dim V = %o at level p (norm 66179) ***\n", Dimension(cur);

// extract a_P and compare to roots {t,t'} (which root, and sign)
if Dimension(cur) gt 0 then
    B := BasisMatrix(cur);
    printf "\n%-7o %-4o %-10o %-12o note\n", "N(P)","p","{a_P}","{t,t'}";
    cnt := 0;
    for tg in targets do
        if Norm(tg[1]+Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        fl, X := IsConsistent(B, B*T); T := 0;
        if not fl then printf "%-7o %-4o V not stable\n", tg[4], tg[5]; continue; end if;
        ev := { e[1] : e in Eigenvalues(X) };
        qr := { r[1] : r in Roots(y^2 - tg[2]*y + tg[3]) };
        printf "%-7o %-4o %-10o %-12o\n", tg[4], tg[5], ev, qr;
        cnt +:= 1; if cnt ge 16 then break; end if;
    end for;
end if;
