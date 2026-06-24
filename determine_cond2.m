// determine_cond2.m -- Decide whether Curve A has bad reduction at the prime above 2.
// This fixes the 2-part of the conductor, hence the level of the 2-dim mod-11 rep W.
// 2 is inert in Q(sqrt(5)), so (2) is a single prime of norm 4.
//
// Run: magma determine_cond2.m

AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, ConductorExponentAt;

k<s5> := QuadraticField(5);
OK := Integers(k);
p2 := Factorization(2*OK)[1][1];
printf "prime above 2: norm %o, ramification %o\n", Norm(p2), RamificationIndex(p2);

QI := [k|
    1,
    1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895)
];
C := Genus2CurveFromIgusa(QI, k);
Cmin, d := MinimalTwist(C, k);
printf "minimal twist d = %o\n", d;

// Integral model and its discriminant valuation at p2.
f0, h0 := HyperellipticPolynomials(Cmin);
m := LCM([Integers()| Denominator(c) : c in Coefficients(f0) cat Coefficients(h0) cat [k!1]]);
R := Parent(f0);
CI := HyperellipticCurve(R!(m^2*f0), R!(m*h0));
disc := Discriminant(CI);
printf "v_{p2}(disc of integral model) = %o\n", Valuation(OK!disc, p2);
printf "v_{p2}(disc) at the *2-part* of Norm: 2-adic val of Norm(disc) = %o\n",
    Valuation(Integers()!Norm(OK!disc), 2);

// Try Magma's conductor exponent at 2 (expected to fail per HANDOFF).
printf "\nConductorExponentAt(Cmin, p2) = ";
e2 := ConductorExponentAt(Cmin, p2);
printf "%o  (-1 means Magma could not compute it)\n", e2;

// Also: does Cmin have good reduction at p2 by the Igusa/discriminant test?
printf "\nIgusa J10 valuation at p2 (J10=0 => bad): v_{p2}(J10) = %o\n",
    Valuation(QI[5], p2);
printf "(J10 is the genus-2 discriminant analogue; v>0 at p2 signals bad reduction)\n";
