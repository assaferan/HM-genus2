// match_twisted.m -- match rho_i for the PRIMITIVE model A_prim = A_recon (x) chi_g (good at
// 499,13711,88301,231481,huge). Fingerprint of rho^prim = chi_g-twist of trnm: tTr -> chi_g(P)*tTr.
// Search untwisted V+ at Pbad*(2)^{0,1,2}.
//
// Run: magma match_twisted.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf, BadPrimesFromInvariants;
Qx<X> := PolynomialRing(Rationals());
K<t> := NumberField(X^2 - X - 1); OK := Integers(K); s5 := 2*t - 1;
PK<x> := PolynomialRing(K);
F11 := GF(11);

LOGF := "match_twisted_out.txt";
PrintFile(LOGF, "# match_twisted" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K); Cint := IntegralModelOf(Cmin);
f, h := HyperellipticPolynomials(Cint);
F0 := PK![OK!c : c in Coefficients(f + (Parent(f)!h)^2/4)];
cont := &+[ideal<OK|c> : c in Coefficients(F0) | c ne 0];
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(cont)]);
_, alpha := IsPrincipal(Msq);
F1 := PK![ (OK!c)/alpha^2 : c in Coefficients(F0) ];
_, g := IsPrincipal(&+[ideal<OK|c> : c in Coefficients(F1) | c ne 0]);
LOG(Sprintf("g (twist) norm = %o", Norm(g)));

Pbad := [P : P in BadPrimesFromInvariants(QI, K) | Norm(P) eq 66179][1];
p2 := Factorization(2*OK)[1][1];

// mod11_trnm was computed over QuadraticField(5) with rid = red(s5). Rebuild primes here.
load "mod11_trnm.txt";
function ChiG(P)
    Fq, red := ResidueClassField(P);
    v := red(OK!g); if v eq 0 then return 0; end if;
    return IsSquare(v) select F11!1 else F11!(-1);
end function;
targets := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    p := row[2]; rid := row[3];
    if rid eq -1 then P := ideal<OK|p>; else P := ideal<OK|p, 2*t-1-rid>; end if;
    cg := ChiG(P);
    if cg eq 0 then continue; end if;
    Append(~targets, <P, cg*F11!row[4], F11!row[5], row[1]>);   // twisted tTr, same tNm
end for;
Sort(~targets, func<a,b | a[4]-b[4]>);

procedure TryLevel(Level, NFIND)
    LOG(Sprintf("== level norm %o ==", Norm(Level)));
    M := HilbertCuspForms(K, Level, [2,2]); d := Dimension(M);
    LOG(Sprintf("  dim %o", d)); if d eq 0 then return; end if;
    Id := IdentityMatrix(F11, d); Vp := VectorSpace(F11, d); Vt := Vp; nf := 0;
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11); T2 := T*T;
        kp := Kernel(T2 - tg[2]*T + tg[3]*Id); km := Kernel(T2 + tg[2]*T + tg[3]*Id);
        Vp := Vp meet kp; Vt := Vt meet (kp+km); T:=0;T2:=0;kp:=0;km:=0; nf+:=1;
        LOG(Sprintf("    N=%o: V+=%o Vtw=%o", tg[4], Dimension(Vp), Dimension(Vt)));
        if Dimension(Vt) eq 0 then break; end if;
        if nf ge NFIND then break; end if;
    end for;
    LOG(Sprintf("  => V+ = %o, Vtw = %o", Dimension(Vp), Dimension(Vt)));
    if Dimension(Vp) gt 0 then LOG("  *** UNTWISTED MATCH (twisted fingerprint)! ***"); end if;
end procedure;

TryLevel(Pbad, 14);
TryLevel(Pbad*p2, 14);
TryLevel(Pbad*p2^2, 12);
LOG("DONE.");
exit;
