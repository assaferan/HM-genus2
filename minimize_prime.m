// minimize_prime.m -- per-prime genus-2 minimization at the 4 residual-content primes
// of Curve A's model (after square-content removal): 499, 13711, 88301, 231481.
// Greedy search over GL2 generators (scale by pi, invert, translate) minimizing v_P(content),
// which for a curve with good reduction at P should reach 0.
//
// Run: magma minimize_prime.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
Qx<X> := PolynomialRing(Rationals());
K<t> := NumberField(X^2 - X - 1); OK := Integers(K); s5 := 2*t - 1;
PK<x> := PolynomialRing(K);

LOGF := "minimize_prime_out.txt";
PrintFile(LOGF, "# minimize_prime" : Overwrite := true);
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

// content valuation at prime P (of the coefficient ideal)
function contVal(F, P)
    return Min([Valuation(OK!c, P) : c in Coefficients(F) | c ne 0]);
end function;
// remove even content globally (square part) to renormalize after a transform
function RenormSquare(F)
    cs := [OK!c : c in Coefficients(F) | c ne 0];
    if #cs eq 0 then return F; end if;
    cn := &+[ideal<OK|c> : c in cs];
    M := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(cn)]);
    ok, a := IsPrincipal(M);
    if ok and a ne 0 then return PK![ (OK!c)/a^2 : c in Coefficients(F) ]; end if;
    return F;
end function;
// apply x -> (a x + b)/(c x + d) to binary sextic F
function Tr(F, a,b,c,d)
    cs := Coefficients(F); while #cs lt 7 do Append(~cs, K!0); end while;
    num := a*x + b; den := c*x + d;
    return &+[ cs[i+1] * num^i * den^(6-i) : i in [0..6] ];
end function;

primes := [499, 13711, 88301, 231481];
for pr in primes do
    Pfac := Factorization(pr*OK);
    for pe in Pfac do
        P := pe[1];
        v0 := contVal(F1, P);
        if v0 eq 0 then continue; end if;   // clean at this P already
        LOG(Sprintf("=== p=%o, prime norm %o, initial content val = %o ===", pr, Norm(P), v0));
        LOG(Sprintf("    Newton (v_P of coeffs c0..c6): %o",
            [Valuation(OK!c, P) : c in ([Coefficient(F1,i) : i in [0..6]])]));
        pi := UniformizingElement(P);
        gens := [ <pi,K!0,K!0,K!1>, <K!1,K!0,K!0,pi>, <K!0,K!1,K!1,K!0>,
                  <K!1,K!1,K!0,K!1>, <K!1,-1,K!0,K!1>, <K!1,t,K!0,K!1>, <K!1,-t,K!0,K!1> ];
        cur := F1; curv := v0; improved := true; steps := 0;
        while improved and steps lt 40 do
            improved := false; steps +:= 1;
            for g in gens do
                G := RenormSquare(Tr(cur, g[1],g[2],g[3],g[4]));
                if &and[ c in OK : c in Coefficients(G) ] and Degree(G) ge 5 then
                    vv := contVal(G, P);
                    if vv lt curv then cur := G; curv := vv; improved := true; break; end if;
                end if;
            end for;
        end while;
        LOG(Sprintf("    after %o steps: content val at P = %o %o",
            steps, curv, curv eq 0 select "(CLEAN)" else "(still nonzero)"));
    end for;
end for;
LOG("DONE.");
exit;
