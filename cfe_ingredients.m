// cfe_ingredients.m -- verify the CFE ingredients for A_true = A_recon (x) chi_g (twisted model,
// good at 499 etc.): (1) is A_true good at 499,13711,88301,231481? (2) the 66179 mult factor.
//
// Run: magma cfe_ingredients.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K);

LOGF := "cfe_ingredients_out.txt";
PrintFile(LOGF, "# cfe_ingredients" : Overwrite := true);
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
_, g := IsPrincipal(&+[ideal<OK|c> : c in Coefficients(F1)|c ne 0]);
Fprim := PK![ c/g : c in Coefficients(F1) ];   // A_true model (twist by g)
LOG(Sprintf("g norm %o; Fprim built (A_true, good at content primes)", Norm(g)));

// per-prime reducer (from mod11_rep): flip+shift to restore degree 6
function TryReduceModP(Fpoly, P)
    Fp, redP := ResidueClassField(P); Rx<xx> := PolynomialRing(Fp);
    n := 6;
    // clear P-denominators: scale Fpoly by pi^k so all coeffs P-integral, min val 0
    vals := [Valuation(c, P) : c in Coefficients(Fpoly) | c ne 0];
    m := Min(vals); pi := UniformizingElement(P);
    Fs := PK![ c/pi^m : c in Coefficients(Fpoly) ];   // note: /pi^m may be a twist if m odd; ok for good-reduction test
    for flip in [0,1] do
        fpoly := flip eq 0 select Fs else PK![Coefficient(Fs, n-i) : i in [0..n]];
        for sh in [0..20] do
            fp := Rx![redP(OK!c) : c in Coefficients(PK!Evaluate(fpoly, x + (OK!sh)))];
            if Degree(fp) eq n and IsSquarefree(fp) then return true, HyperellipticCurve(fp); end if;
        end for;
    end for;
    return false, _;
end function;

for pr in [7, 13, 499, 13711, 88301, 231481] do
    for pe in Factorization(pr*OK) do
        P := pe[1];
        ok, Cbar := TryReduceModP(Fprim, P);
        if ok then
            LOG(Sprintf("  p=%o (norm %o): A_true GOOD, L=%o", pr, Norm(P), LPolynomial(Cbar)));
        else
            LOG(Sprintf("  p=%o (norm %o): A_true reduction NOT smooth (bad or hard)", pr, Norm(P)));
        end if;
    end for;
end for;

// 66179 factor (mult reduction): reduce Fprim mod P66179, nodal -> elliptic normalization
P66 := [pe[1] : pe in Factorization(66179*OK) | Norm(pe[1]) eq 66179][1];
Fp, red := ResidueClassField(P66); Rq<xx> := PolynomialRing(Fp);
vals := [Valuation(c,P66) : c in Coefficients(Fprim)|c ne 0]; m := Min(vals); pi:=UniformizingElement(P66);
fbar := Rq![red(OK!(c/pi^m)) : c in Coefficients(Fprim)];
LOG(Sprintf("66179: fbar deg %o, sqfree %o", Degree(fbar), IsSquarefree(fbar)));
gg := GCD(fbar, Derivative(fbar));
LOG(Sprintf("  gcd(fbar,fbar') deg = %o (1 => single node/mult red)", Degree(gg)));
LOG("DONE.");
exit;
