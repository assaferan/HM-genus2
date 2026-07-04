// cfe_final.m -- determine c = a_{p2}(A_true) via the functional equation.
// A_true = A_recon (x) chi_g, conductor p * (2)^c, so N(W)=5^4*66179*4^c. Search c.
// Euler factors via ZoomL (direct reduction + zoom for non-minimal primes), twisted by chi_g.

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K);
ZT<T> := PolynomialRing(Integers());
LOGF := "cfe_final_out.txt"; PrintFile(LOGF, "# cfe_final" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787), 1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089), 1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K); Cint := IntegralModelOf(Cmin);
ff,hh := HyperellipticPolynomials(Cint);
F0 := PK![OK!c : c in Coefficients(ff + (Parent(ff)!hh)^2/4)];
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(&+[ideal<OK|c> : c in Coefficients(F0)|c ne 0])]);
_, alpha := IsPrincipal(Msq);
F1 := PK![ (OK!c)/alpha^2 : c in Coefficients(F0) ];
_, g := IsPrincipal(&+[ideal<OK|c> : c in Coefficients(F1)|c ne 0]);
Fprim := PK![ c/g : c in Coefficients(F1) ];
LOG("setup done");

function ChiG(P) Fp,red := ResidueClassField(P); v := red(OK!g); if v eq 0 then return 0; end if; return IsSquare(v) select 1 else -1; end function;

// reduce model F at P (zoom for clustered/non-minimal); return LPolynomial in T (degree<=4), or 0
function ZoomL(F, P)
    Fp, red := ResidueClassField(P); Rx<u> := PolynomialRing(Fp);
    // FAST path: direct reduction (works when F is P-integral, e.g. F1)
    isint := true; fbar := Rx!0;
    try fbar := Rx![ red(OK!c) : c in Coefficients(F) ]; catch e; isint := false; end try;
    if isint and fbar ne 0 and Degree(fbar) in {5,6} and IsSquarefree(fbar) then
        return ZT!LPolynomial(HyperellipticCurve(fbar));
    end if;
    // SLOW path: content-clearing + zoom
    _, pi := IsPrincipal(P);
    for depth in [0..15] do
        cs := [c : c in Coefficients(F) | c ne 0]; if #cs eq 0 then return 0; end if;
        m := Min([Valuation(c,P) : c in cs]);
        if m ne 0 then F := PK![ c/pi^(2*(m div 2)) : c in Coefficients(F) ]; end if;
        fbar := Rx![ red(OK!c) : c in Coefficients(F) ];
        if fbar eq 0 then return 0; end if;
        if Degree(fbar) in {5,6} and IsSquarefree(fbar) then return ZT!LPolynomial(HyperellipticCurve(fbar)); end if;
        gd := GCD(fbar, Derivative(fbar));
        if Degree(gd) eq 0 then F := PK![Coefficient(F,6-i) : i in [0..6]]; continue; end if;
        rr := Roots(gd); if #rr eq 0 then return 0; end if;
        r := OK!(rr[1][1] @@ red);
        F := Evaluate(F, pi*x + r);
    end for;
    return 0;
end function;

CUTOFF := 47000;
Ltrue := AssociativeArray(); nfail := 0;
for p in PrimesUpTo(CUTOFF) do
    if p eq 2 then continue; end if;
    for pe in Factorization(p*OK) do
        P := pe[1]; fP := Degree(ResidueClassField(P));
        if fP eq 1 then _r,_d := ResidueClassField(P); rid := Integers()!_d(OK!s5); else rid := -1; end if;
        cg := ChiG(P);
        L := ZoomL(F1, P);
        if L eq 0 then
            L := ZoomL(Fprim, P);                    // A_true directly (recon-bad primes)
            if L ne 0 then Ltrue[<p,rid>] := L; continue; end if;   // already true, no twist
            nfail +:= 1; LOG(Sprintf("  FAIL p=%o rid=%o norm=%o", p, rid, Norm(P))); continue;
        end if;
        Ltrue[<p,rid>] := Evaluate(L, cg*T);         // A_recon L, twisted to A_true
    end for;
end for;
LOG(Sprintf("precompute: %o factors, %o failures", #Keys(Ltrue), nfail));

function LocalFactor(p, deg)
    if p eq 2 then return ZT!1; end if;
    res := ZT!1;
    for pe in Factorization(p*OK) do
        P := pe[1]; fP := Degree(ResidueClassField(P));
        if fP eq 1 then _r,_d := ResidueClassField(P); rid := Integers()!_d(OK!s5); else rid := -1; end if;
        if not IsDefined(Ltrue, <p,rid>) then return ZT!1; end if;
        res := res * Evaluate(Ltrue[<p,rid>], T^fP);
    end for;
    return res mod T^(deg+1);
end function;
cf := func< p, d | LocalFactor(p, d) >;

LOG("CFE search:");
for c in [5..12] do
    N := 5^4 * 66179 * 4^c;
    L := LSeries(2, [0,0,0,0,1,1,1,1], N, cf : Sign := 0, Precision := 1);
    try err := CheckFunctionalEquation(L); LOG(Sprintf("  c=%o (N=%o): |CFE|=%o", c, N, err));
    catch e; LOG(Sprintf("  c=%o: err %o", c, e`Object)); end try;
end for;
LOG("DONE."); exit;
