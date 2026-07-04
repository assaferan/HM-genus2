// minimize_hard.m -- get A_recon's good Euler factor at the 4 badly-non-minimal primes
// (3 inert, 5 ramified, one prime|19, one prime|191) via greedy v_P(disc) minimization.
//
// Run: magma minimize_hard.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K);

LOGF := "minimize_hard_out.txt";
PrintFile(LOGF, "# minimize_hard" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K); Cint := IntegralModelOf(Cmin);
f, h := HyperellipticPolynomials(Cint);
F0 := PK![OK!c : c in Coefficients(f + (Parent(f)!h)^2/4)];
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(&+[ideal<OK|c> : c in Coefficients(F0)|c ne 0])]);
_, alpha := IsPrincipal(Msq);
F1 := PK![ (OK!c)/alpha^2 : c in Coefficients(F0) ];

function Tr(F, a,b,c,d)
    cs := Coefficients(F); while #cs lt 7 do Append(~cs, K!0); end while;
    num := a*x+b; den := c*x+d;
    return &+[ cs[i+1]*num^i*den^(6-i) : i in [0..6] ];
end function;
function RenormSq(F, P)
    cs := [c : c in Coefficients(F) | c ne 0]; if #cs eq 0 then return F; end if;
    vv := [Valuation(c, P) : c in cs]; m := Min(vv);
    if m eq 0 then return F; end if;
    pi := UniformizingElement(P); e := 2*(m div 2);   // remove even part only (isomorphism)
    return PK![ c/pi^e : c in Coefficients(F) ];
end function;
function vDisc(F, P)
    d := Discriminant(F); if d eq 0 then return 10^9; end if; return Valuation(d, P);
end function;

function MinimizeAt(F, P)
    pi := UniformizingElement(P);
    reps := [OK| 0,1,-1,2,-2,s5,-s5,1+s5,1-s5,3,-3,s5+2,s5-2];
    gens := [<pi,K!0,K!0,K!1>, <K!1,K!0,K!0,pi>, <K!0,K!1,K!1,K!0>]
            cat [<K!1, K!r, K!0, K!1> : r in reps];
    cur := RenormSq(F, P); bestv := vDisc(cur, P); improved := true; steps := 0;
    while improved and steps lt 60 do
        improved := false; steps +:= 1;
        for gg in gens do
            G := RenormSq(Tr(cur, gg[1],gg[2],gg[3],gg[4]), P);
            if Degree(G) lt 5 then continue; end if;
            v := vDisc(G, P);
            if v lt bestv then cur := G; bestv := v; improved := true; break; end if;
        end for;
    end while;
    return cur, bestv;
end function;

targets := [<3,-1>, <5,0>, <19,10>, <191,177>];
for tp in targets do
    p := tp[1]; rid := tp[2];
    P := (rid eq -1) select ideal<OK|p> else ideal<OK|p, s5-rid>;
    Gmin, v := MinimizeAt(F1, P);
    LOG(Sprintf("p=%o rid=%o: min v_P(disc) = %o", p, rid, v));
    Fp, red := ResidueClassField(P); Rx<xx> := PolynomialRing(Fp);
    // clear content and reduce
    vv := [Valuation(c,P) : c in Coefficients(Gmin)|c ne 0]; m := Min(vv); pi := UniformizingElement(P);
    fbar := Rx![red(OK!(c/pi^m)) : c in Coefficients(Gmin)];
    if Degree(fbar) eq 6 and IsSquarefree(fbar) then
        LOG(Sprintf("  GOOD reduction; L_recon = %o", LPolynomial(HyperellipticCurve(fbar))));
    elif Degree(fbar) eq 5 and IsSquarefree(fbar) then
        LOG(Sprintf("  GOOD (deg 5); L_recon = %o", LPolynomial(HyperellipticCurve(fbar))));
    else
        LOG(Sprintf("  still not smooth: deg %o sqfree %o", Degree(fbar), IsSquarefree(fbar)));
    end if;
end for;
LOG("DONE.");
exit;
