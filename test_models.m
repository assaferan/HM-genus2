// test_models.m -- do F1 or f_deg6 reduce smoothly at the 4 hard primes 3,5,19,191?
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K);
LOGF := "test_models_out.txt"; PrintFile(LOGF, "# test_models" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;
QI := [K| 1, 1/1459240*(243125*s5 - 482787), 1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089), 1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K);
f0 := HyperellipticPolynomials(Cmin);
Rk<xx> := PolynomialRing(K);
function Mob(f,a,b,c,d) n:=Degree(f); return &+[Coefficient(f,i)*(a*xx+b)^i*(c*xx+d)^(n-i) : i in [0..n]]; end function;
f_mob := Mob(f0, 1, -s5+2, K!0, s5+2);
f_deg6 := Mob(f_mob, -5*s5-5, -5*s5+3, K!0, (-5*s5-11)/2);
// also F1
Cint := IntegralModelOf(Cmin); ff,hh := HyperellipticPolynomials(Cint);
F0 := PK![OK!c : c in Coefficients(ff + (Parent(ff)!hh)^2/4)];
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(&+[ideal<OK|c> : c in Coefficients(F0)|c ne 0])]);
_, alpha := IsPrincipal(Msq);
F1 := PK![ (OK!c)/alpha^2 : c in Coefficients(F0) ];

function TryModel(F, P)   // reduce model F at P (even-content only), try shifts; return L or 0
    Fp, red := ResidueClassField(P); Rx<u> := PolynomialRing(Fp);
    cs := [c : c in Coefficients(F) | c ne 0]; if #cs eq 0 then return 0; end if;
    m := Min([Valuation(c,P) : c in cs]);
    if IsOdd(m) then return 0; end if;   // odd content would twist
    pi := UniformizingElement(P);
    fbar0 := Rx![red(OK!(c/pi^m)) : c in Coefficients(F)];
    lifts := [Fp| ]; for a in Fp do Append(~lifts, a); end for;
    for flip in [0,1] do
        base := flip eq 0 select fbar0 else Rx![Coefficient(fbar0, 6-i) : i in [0..6]];
        for s in lifts do
            g := Evaluate(base, u + s);
            if Degree(g) in {5,6} and IsSquarefree(g) then return LPolynomial(HyperellipticCurve(g)); end if;
        end for;
    end for;
    return 0;
end function;

for tp in [<3,-1>,<5,0>,<19,10>,<191,177>] do
    p := tp[1]; rid := tp[2];
    P := (rid eq -1) select ideal<OK|p> else ideal<OK|p, s5-rid>;
    L1 := TryModel(F1, P); L2 := TryModel(f_deg6, P);
    LOG(Sprintf("p=%o: F1 -> %o ;  f_deg6 -> %o", p, L1 eq 0 select "FAIL" else Sprint(L1), L2 eq 0 select "FAIL" else Sprint(L2)));
end for;
LOG("DONE."); exit;
