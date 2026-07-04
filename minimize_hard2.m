// minimize_hard2.m -- resolve badly-non-minimal primes by zoom-search: roots cluster P-adically;
// x -> pi*x + r (r a repeated root mod P, pi = P-local generator) separates them. A_recon's L.
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K);
LOGF := "minimize_hard2_out.txt"; PrintFile(LOGF, "# minimize_hard2" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;
QI := [K| 1, 1/1459240*(243125*s5 - 482787), 1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089), 1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K); Cint := IntegralModelOf(Cmin);
ff,hh := HyperellipticPolynomials(Cint);
F0 := PK![OK!c : c in Coefficients(ff + (Parent(ff)!hh)^2/4)];
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(&+[ideal<OK|c> : c in Coefficients(F0)|c ne 0])]);
_, alpha := IsPrincipal(Msq);
F1 := PK![ (OK!c)/alpha^2 : c in Coefficients(F0) ];

function ZoomL(F, P)
    _, pi := IsPrincipal(P);                        // P-LOCAL uniformiser (unit at all other primes)
    Fp, red := ResidueClassField(P); Rx<u> := PolynomialRing(Fp);
    for depth in [0..6] do
        // clear even content, reduce
        cs := [c : c in Coefficients(F) | c ne 0]; if #cs eq 0 then return 0; end if;
        m := Min([Valuation(c,P) : c in cs]);
        F := PK![ c / pi^(2*(m div 2)) : c in Coefficients(F) ];
        fbar := Rx! [ red(OK!c) : c in Coefficients(F) ];
        if Degree(fbar) in {5,6} and IsSquarefree(fbar) then return LPolynomial(HyperellipticCurve(fbar)); end if;
        gd := GCD(fbar, Derivative(fbar));
        if Degree(gd) eq 0 then
            F := PK![Coefficient(F,6-i) : i in [0..6]];   // repeated root at infinity: flip
            continue;
        end if;
        rr := Roots(gd); r := OK ! (rr[1][1] @@ red);
        F := Evaluate(F, pi*x + r);                       // zoom into the cluster at r
    end for;
    return 0;
end function;

for tp in [<3,-1>,<5,0>,<19,10>,<191,177>] do
    p := tp[1]; rid := tp[2];
    P := (rid eq -1) select ideal<OK|p> else ideal<OK|p, s5-rid>;
    L := ZoomL(F1, P);
    LOG(Sprintf("p=%o rid=%o: %o", p, rid, L eq 0 select "FAIL" else Sprint(L)));
end for;
LOG("DONE."); exit;
