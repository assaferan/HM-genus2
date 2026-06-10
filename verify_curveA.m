// verify_curveA.m -- Recreate and verify Curve A from the Q(sqrt(5)) (1,11) search.
//
// Steps:
//   1. Paste the Igusa invariants.
//   2. Call Genus2CurveFromIgusa to get a Weierstrass model.
//   3. Apply MinimalTwist (should give d=1 for Curve A).
//   4. Verify: bad prime N(p)=66179 has conductor exponent 1,
//              both p|11 have conductor exponent 0.
//   5. Spot-check a few Frobenius L-polynomials against data_qsqrt5_gp11.txt.
//
// Run: magma verify_curveA.m

AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, ConductorExponentAt, BadPrimesFromInvariants, PrimesAbove;

k<s5> := QuadraticField(5);
OK := Integers(k);
ZT<T> := PolynomialRing(Integers());

// -----------------------------------------------------------------------
// Igusa invariants for Curve A (= survivors 1 and 2 in survivors_qsqrt5_gp11.m)
// -----------------------------------------------------------------------
QI := [k|
    1,
    1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895)
];
printf "Curve A over K = Q(sqrt(5)):\n";
printf "  J4  = %o\n  J10 = %o\n\n", QI[2], QI[5];

// -----------------------------------------------------------------------
// Reconstruct Weierstrass model
// -----------------------------------------------------------------------
printf "Step 1: Reconstruct curve from Igusa invariants...\n";
C := Genus2CurveFromIgusa(QI, k);
f, h := HyperellipticPolynomials(C);
// Show only the degree and leading coefficient (coefficients are 100+ digits)
printf "  Weierstrass model: y^2 = f(x)  [h=%o]\n", h;
printf "  Degree of f: %o,  leading coeff has norm ~ 10^%o\n",
    Degree(f), #IntegerToString(Ceiling(Abs(Norm(LeadingCoefficient(f)))));
printf "\n";

// -----------------------------------------------------------------------
// Minimal twist
// -----------------------------------------------------------------------
printf "Step 2: Minimal quadratic twist...\n";
Cmin, d := MinimalTwist(C, k);
printf "  Twist factor d = %o  (expect 1 for Curve A)\n\n", d;

// -----------------------------------------------------------------------
// Bad primes from Igusa support
// -----------------------------------------------------------------------
printf "Step 3: Bad primes (Igusa valuation criterion, exact):\n";
bad := BadPrimesFromInvariants(QI, k);
for p in bad do
    printf "  N(p) = %o  (rational prime %o)\n", Norm(p), Minimum(p);
end for;
printf "\n";

// -----------------------------------------------------------------------
// Conductor exponents
// -----------------------------------------------------------------------
printf "Step 4: Conductor exponents:\n";
printf "  At the bad prime:\n";
for p in bad do
    e := ConductorExponentAt(Cmin, p);
    printf "    N(p) = %o: exponent = %o  [%o]  (expect 1)\n",
        Norm(p), e, e eq 1 select "OK" else "UNEXPECTED";
end for;
printf "  At primes above 11 (11 splits in K):\n";
for p11 in PrimesAbove(k, 11) do
    e := ConductorExponentAt(Cmin, p11);
    printf "    N(p11) = %o: exponent = %o  [%o]  (expect 0)\n",
        Norm(p11), e, e eq 0 select "OK" else "UNEXPECTED";
end for;
printf "\n";

// -----------------------------------------------------------------------
// Frobenius spot-check (using the same TryReduceModP as data_qsqrt5_gp11.m)
// -----------------------------------------------------------------------
function MakeIntegral(C)
    f0, h0 := HyperellipticPolynomials(C);
    K0 := BaseRing(C); OK0 := Integers(K0);
    m := LCM([Integers()| Denominator(c) : c in Coefficients(f0) cat Coefficients(h0) cat [K0!1]]);
    return HyperellipticCurve(Parent(f0)!(m^2*f0), Parent(f0)!(m*h0));
end function;

// Exactly the version from data_qsqrt5_gp11.m (integer shifts 0..20 + flip).
function TryReduceModP(CI, P)
    K0 := BaseRing(CI); OK0 := Integers(K0);
    Fp, redP := ResidueClassField(P); Rx<x> := PolynomialRing(Fp);
    f0, _ := HyperellipticPolynomials(CI); n := Degree(f0); R := Parent(f0);
    for flip in [0, 1] do
        f := flip eq 0 select f0
             else Polynomial([Coefficient(f0, n-i) : i in [0..n]]);
        for sh in [0..20] do
            fp := Rx![redP(OK0!c) : c in Coefficients(Evaluate(f, R.1 + sh))];
            if Degree(fp) eq n then return true, HyperellipticCurve(fp, Rx!0); end if;
        end for;
    end for;
    return false, _;
end function;

CI := MakeIntegral(Cmin);
bad_set := {p : p in bad};
R11<t11> := PolynomialRing(GF(11));

printf "Step 5: Frobenius L-polynomials (spot-check vs data_qsqrt5_gp11.txt):\n";
printf "  Expected at N(p)=49 (p=7, inert): 2401*T^4 + 98*T^3 - 58*T^2 + 2*T + 1\n";
printf "  Expected at N(p)=11 (first):       121*x^4 + 33*x^3 + 18*x^2 + 3*x + 1\n\n";
printf "  %-6o  %-5o  %o\n", "N(p)", "p", "L_p(T)";
printf "  " cat "-"^75 cat "\n";

count := 0;
pp := 3;
while count lt 10 and pp lt 200 do
    for pid in PrimesAbove(k, pp) do
        if pid in bad_set then continue; end if;
        if ConductorExponentAt(Cmin, pid) ne 0 then continue; end if;
        ok, Cp := TryReduceModP(CI, pid);
        if not ok then
            printf "  %-6o  %-5o  [reduction failed]\n", Norm(pid), pp;
        else
            lp := ZT ! LPolynomial(Cp);
            fac := Factorization(R11 ! lp);
            printf "  %-6o  %-5o  %o\n", Norm(pid), pp, lp;
            count +:= 1;
        end if;
        if count ge 10 then break; end if;
    end for;
    pp := NextPrime(pp);
end while;

printf "\n=== Verification summary ===\n";
printf "Bad prime N(p)=66179: conductor exponent %o  (expected 1)\n",
    ConductorExponentAt(Cmin, bad[1]);
printf "Both primes above 11: exponents %o, %o  (both expected 0)\n",
    ConductorExponentAt(Cmin, PrimesAbove(k,11)[1]),
    ConductorExponentAt(Cmin, PrimesAbove(k,11)[2]);
