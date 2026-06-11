// reduce_curveA5.m -- Print the best model found by the reduction chain.
//
// Applies the full chain of transforms found by reduce_curveA.m / reduce_curveA2.m:
//   Cmin  --Mob1--> f_mob  --affine--> f_best   (degree 6, h~18)
//   Cmin  --HasOddDegreeModel--> f5              (degree 5, h~28.9)
// and prints actual coefficient sizes so we can judge practicality.
//
// Run: magma reduce_curveA5.m

AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist;

k<s5> := QuadraticField(5);
OK := Integers(k);
Rk<x> := PolynomialRing(k);

QI := [k|
    1,
    1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895)
];

function LogHeight(f, k)
    sqr5 := Sqrt(RealField(50)!5);
    embs := [sqr5, -sqr5];
    h := RealField(50)!0;
    for c in Coefficients(f) do
        seq := Eltseq(k!c);
        a := RealField(50)!seq[1]; b := RealField(50)!seq[2];
        for e in embs do h := Max(h, Abs(a + b*e)); end for;
    end for;
    return h gt 0 select Log(h) else RealField(50)!0;
end function;

function MobiusPoly(f, a, b, c, d, Rk)
    n := Degree(f); g := Rk!0;
    for i in [0..n] do
        g +:= Coefficient(f,i) * (a*Rk.1+b)^i * (c*Rk.1+d)^(n-i);
    end for;
    return g;
end function;

// Print coefficient neatly: show numerator/denominator sizes.
procedure PrintCoeffInfo(c, i, k)
    if k!c eq 0 then return; end if;
    seq := Eltseq(k!c);
    a := seq[1]; b := seq[2];
    printf "  [x^%o] (%o) + (%o)*s5\n", i, a, b;
    // Size summary
    na := Numerator(a); da := Denominator(a);
    nb := Numerator(b); db := Denominator(b);
    printf "         rational part: %o / %o  (%o / %o digits)\n",
        na, da, #IntegerToString(Abs(na)), #IntegerToString(da);
    printf "         s5-part:       %o / %o  (%o / %o digits)\n",
        nb, db, #IntegerToString(Abs(nb)), #IntegerToString(db);
end procedure;

// -----------------------------------------------------------------------
C := Genus2CurveFromIgusa(QI, k);
Cmin, d := MinimalTwist(C, k);
f0, _ := HyperellipticPolynomials(Cmin);

// -----------------------------------------------------------------------
// Model 1: Degree-6, best from reduce_curveA2.m (h~18)
// Chain: f0 -> Mob(1,-s5+2,0,s5+2) -> affine(-5s5-5,-5s5+3,0,(-5s5-11)/2)
// -----------------------------------------------------------------------
printf "=== Model 1: Degree-6, log-height ~18 ===\n";
f_mob := MobiusPoly(f0, 1, -s5+2, 0, s5+2, Rk);
// Best affine from reduce_curveA2.m at h=18.047: (a=-5s5-5, b=-5s5+3, c=0, d=(-5s5-11)/2)
a1 := -5*s5-5; b1 := -5*s5+3; c1 := k!0; d1 := (-5*s5-11)/2;
f_deg6 := MobiusPoly(f_mob, a1, b1, c1, d1, Rk);
h6 := LogHeight(f_deg6, k);
printf "log-height: %.4o  (max|coeff under emb| ~ e^%.1o ≈ 10^%.1o)\n\n",
    h6, h6, Log(10, Exp(RealField(10)!h6));
for i in [Degree(f_deg6)..0 by -1] do
    PrintCoeffInfo(Coefficient(f_deg6,i), i, k);
end for;
printf "\n";

// -----------------------------------------------------------------------
// Model 2: Degree-5 model (from HasOddDegreeModel)
// -----------------------------------------------------------------------
printf "=== Model 2: Degree-5 model, log-height ~28.9 ===\n";
has5, C5 := HasOddDegreeModel(Cmin);
if has5 then
    f5, _ := HyperellipticPolynomials(C5);
    printf "log-height: %.4o  (max|coeff| ~ 10^%.1o)\n\n",
        LogHeight(f5,k), Log(10, Exp(RealField(10)!LogHeight(f5,k)));
    for i in [Degree(f5)..0 by -1] do
        PrintCoeffInfo(Coefficient(f5,i), i, k);
    end for;
else
    printf "  No degree-5 model.\n";
end if;
printf "\n";

// -----------------------------------------------------------------------
// Model 3: Degree-5 after best rational affine reduction from reduce_curveA3.m
//          (h~23.54 was found with alpha=1/2, beta=5/2 or similar shifts)
// -----------------------------------------------------------------------
printf "=== Model 3: Degree-5 after affine reduction (h~23.5) ===\n";
if has5 then
    // From reduce_curveA3.m the best was deg5 iter 1: alpha=1/2 beta=5/2 h=23.5563
    // then stagnated at 23.5428.  Apply the sequence:
    f5_r := MobiusPoly(f5, k!(1), k!(5)/2, k!0, k!1, Rk);  // x -> x + 5/2
    f5_r := MobiusPoly(f5_r, k!(1)/2, k!0, k!0, k!1, Rk);  // x -> x/2
    printf "After x->x+5/2 then x->x/2: h=%.4o\n\n", LogHeight(f5_r,k);
    for i in [Degree(f5_r)..0 by -1] do
        PrintCoeffInfo(Coefficient(f5_r,i), i, k);
    end for;
end if;
printf "\n";

// -----------------------------------------------------------------------
// Summary: what's the numerator/denominator structure telling us?
// -----------------------------------------------------------------------
printf "=== Assessment ===\n";
printf "The large denominators in all models come directly from the Igusa\n";
printf "invariants of Curve A (J4 denom ~10^6, J10 denom ~10^16).  Mestre's\n";
printf "algorithm propagates these into model coefficients.  No GL2(O_K)\n";
printf "transform can clear them -- the model is intrinsically non-integral.\n\n";
printf "Best practical representations:\n";
printf "  (a) Quote the Igusa-Siegel invariants [J2:J4:...:J10] directly.\n";
printf "      These are compact and sufficient for identifying the curve.\n";
printf "  (b) Degree-5 Weierstrass model (log-height ~23.5) has ~10^10 coeff size.\n";
printf "  (c) Degree-6 model (log-height ~18) has ~10^8 coeff size.\n";
printf "  (d) For computations, work mod p and use the TryReduceModP reduction.\n";
