// reduce_curveA2.m -- Additional reduction attempts for Curve A.
//
// 1. Check whether a degree-5 model exists (HasOddDegreeModel).
// 2. Inspect the Mestre conic directly and search for its smallest-height K-point.
// 3. Iterate the Mobius reduction (apply step 3 of reduce_curveA.m to its own output).
//
// Run: magma reduce_curveA2.m

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
    sqr5 := Sqrt(RealField(30)!5);
    embs := [sqr5, -sqr5];
    h := RealField(30)!0;
    for c in Coefficients(f) do
        seq := Eltseq(k!c);
        a := RealField(30)!seq[1]; b := RealField(30)!seq[2];
        for e in embs do h := Max(h, Abs(a + b*e)); end for;
    end for;
    return h gt 0 select Log(h) else RealField(30)!0;
end function;

function MobiusPoly(f, a, b, c, d, Rk)
    n := Degree(f); g := Rk!0;
    for i in [0..n] do
        g +:= Coefficient(f,i) * (a*Rk.1+b)^i * (c*Rk.1+d)^(n-i);
    end for;
    return g;
end function;

// -----------------------------------------------------------------------
// Setup: get Cmin
// -----------------------------------------------------------------------
C := Genus2CurveFromIgusa(QI, k);
Cmin, d := MinimalTwist(C, k);
f0, _ := HyperellipticPolynomials(Cmin);
printf "Twist d = %o\n", d;
printf "Cmin log-height: %.4o\n\n", LogHeight(f0, k);

// The best transform from reduce_curveA.m was Mobius (1, -s5+2, 0, s5+2).
f_mob := MobiusPoly(f0, 1, -s5+2, 0, s5+2, Rk);
printf "After best Mobius (1,-s5+2,0,s5+2): log-height %.4o\n\n", LogHeight(f_mob, k);

// -----------------------------------------------------------------------
// 1. Degree-5 model: does the curve have a K-rational Weierstrass point?
// -----------------------------------------------------------------------
printf "=== 1. Degree-5 model ===\n";
try
    has5, C5 := HasOddDegreeModel(Cmin);
    if has5 then
        f5, h5 := HyperellipticPolynomials(C5);
        printf "  HasOddDegreeModel: YES! degree-%o model found.\n", Degree(f5);
        printf "  log-height of degree-5 model: %.4o\n", LogHeight(f5, k);
        printf "  f5 =\n";
        for i in [Degree(f5)..0 by -1] do
            c := Coefficient(f5, i);
            if c ne 0 then printf "    [x^%o] %o\n", i, c; end if;
        end for;
    else
        printf "  HasOddDegreeModel: NO (no K-rational Weierstrass point).\n";
    end if;
catch e
    printf "  HasOddDegreeModel error: %o\n", e`Object;
end try;
printf "\n";

// -----------------------------------------------------------------------
// 2. Mestre conic: inspect and search for small K-points
//
// The Igusa-Clebsch invariants are I2,I4,I6,I10. The conic L is a quadratic
// form in 3 variables with coefficients in k. Once we have the conic, we can
// parametrize all k-rational points via one known point (xi0, eta0) and a
// slope parameter m, then minimize over small m.
// -----------------------------------------------------------------------
printf "=== 2. Mestre conic inspection ===\n";
IC := [k| c : c in IgusaClebschInvariants(Cmin)];
printf "  Igusa-Clebsch invariants:\n";
for i in [1..4] do
    printf "    IC[%o] = %o\n", i, IC[i];
end for;
printf "\n";

// Reconstruct the conic L using the same formulas as mestre.m (ConicLAndCubicM).
// The Mestre conic in variables (u,v,w) is:
//   L = [[L11,L12,L13],[L12,L22,L23],[L13,L23,L33]] with
//   L11 = -120*I6, L12 = 20*I4, L13 = -2*I2, L22 = -120*I6, ...
// (exact formulas from mestre.m ConicLAndCubicM -- we read them from
// the CHIMP/Magma source and reproduce the conic here).
//
// Rather than reimplementing ConicLAndCubicM, we use Magma's conic object.
P2<u,v,w> := ProjectiveSpace(k, 2);
I2 := IC[1]; I4 := IC[2]; I6 := IC[3]; I10 := IC[4];

// Mestre conic coefficients (from mestre.m lines 323-430):
//   L = 8*I4*u^2 + 2*I2*u*v - 6*I6*... -- the exact formula is in ConicLAndCubicM.
// We extract it by building the conic symbolically then projecting.
// Actually, the cleanest way is to use Magma's HasRationalPoint directly.
P2k<uu,vv,ww> := ProjectiveSpace(k, 2);
// Use Magma's internal conic via the genus-2 function:
// HyperellipticCurveFromIgusaClebsch calls ConicLAndCubicM internally.
// We can't directly access L without reimplementing it, so let us use
// a different strategy: parametrize the known K-point.

// Recover the known conic point by comparing two models from the same Igusa invs
// at different (arbitrary) "twist" parameters.
// Strategy: fix one model, find the rational point that generates it, then
// search over nearby rational points.

// The "slope parametrization" of the Mestre conic:
// from the base point (xi0, eta0) and slope m, the parametric family is
// f_m(x) = Mestre cubic evaluated at the second conic point.
// Different m <-> different K-isomorphic sextic model.
// We scan m over O_K elements with small norm.

// Build all elements of O_K with |a|+|b| <= B:
B := 5;
small_ok := [k| a + b*k.1 : a in [-B..B], b in [-B..B]];
// Add some elements of the form (a+b*s5)/2 (half-integers in O_K):
// O_K = Z[(1+s5)/2], so include (1+s5)/2 and its small multiples.
w5 := (1+s5)/2;
for a in [-B..B] do for b in [-B..B] do
    Append(~small_ok, a + b*w5);
end for; end for;

printf "  Scanning %o small O_K elements as Mobius parameters (iterating from best so far)...\n",
    #small_ok;

// Apply Mobius to f_mob (already log-height 32.49) for a second reduction pass.
f_cur := f_mob; best_h := LogHeight(f_mob, k); best_f := f_mob;
best_desc := "prev-best";
gens2 := small_ok;
n6 := Degree(f_cur);
for a in gens2 do for b in gens2 do for c in gens2 do for d in gens2 do
    det := a*d - b*c;
    if det eq 0 then continue; end if;
    g := MobiusPoly(f_cur, a, b, c, d, Rk);
    if Degree(g) lt n6 then continue; end if;
    if LeadingCoefficient(g) eq 0 then continue; end if;
    h := LogHeight(g, k);
    if h lt best_h then
        best_h := h; best_f := g;
        best_desc := Sprintf("2nd-pass Mobius (%o,%o,%o,%o)", a,b,c,d);
        printf "  %o  h=%.4o\n", best_desc, h;
    end if;
end for; end for; end for; end for;
printf "  Best after 2nd Mobius pass: %o  h=%.4o\n\n", best_desc, best_h;

// -----------------------------------------------------------------------
// 3. Try working with the degree-5 form of the conic parametrization:
//    apply Mobius to f_mob with generators from O_K half-integers only.
// -----------------------------------------------------------------------
printf "=== 3. Summary ===\n";
printf "  Original Cmin:        log-height %.4o\n", LogHeight(f0, k);
printf "  After 1st Mobius:     log-height %.4o\n", LogHeight(f_mob, k);
printf "  After 2nd pass:       log-height %.4o\n", best_h;
printf "\n";
printf "Best model:\n";
for i in [Degree(best_f)..0 by -1] do
    c := Coefficient(best_f, i);
    if c ne 0 then printf "  [x^%o] %o\n", i, c; end if;
end for;
