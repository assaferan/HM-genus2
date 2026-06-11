// reduce_curveA4.m -- Get the best possible degree-5 model for Curve A.
//
// Strategy: HasOddDegreeModel gives a degree-5 Weierstrass model.
// We then apply the full Mobius search (all generators over a wide O_K set)
// to find the smallest-height degree-5 model.
//
// Run: magma reduce_curveA4.m

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

function MaxCoeffNorm(f, k)
    // max over coefficients of N(c)^{1/2} (geometric mean of embeddings)
    sqr5 := Sqrt(RealField(50)!5);
    best := RealField(50)!0;
    for c in Coefficients(f) do
        seq := Eltseq(k!c);
        a := RealField(50)!seq[1]; b := RealField(50)!seq[2];
        v1 := Abs(a + b*sqr5); v2 := Abs(a - b*sqr5);
        best := Max(best, v1*v2);  // |Norm(c)|
    end for;
    return best;
end function;

function MobiusPoly(f, a, b, c, d, Rk)
    n := Degree(f); g := Rk!0;
    for i in [0..n] do
        g +:= Coefficient(f,i) * (a*Rk.1+b)^i * (c*Rk.1+d)^(n-i);
    end for;
    return g;
end function;

// -----------------------------------------------------------------------
// Get degree-5 model
// -----------------------------------------------------------------------
printf "=== Get degree-5 model ===\n";
C := Genus2CurveFromIgusa(QI, k);
Cmin, d := MinimalTwist(C, k);
printf "Twist d = %o\n", d;

has5, C5 := HasOddDegreeModel(Cmin);
if not has5 then error "No degree-5 model found!"; end if;
f5, _ := HyperellipticPolynomials(C5);
printf "Degree-5 model: y^2 = f(x), deg(f) = %o, log-height = %.4o\n\n",
    Degree(f5), LogHeight(f5, k);

// Print actual coefficients
printf "Raw degree-5 coefficients:\n";
for i in [Degree(f5)..0 by -1] do
    c := Coefficient(f5, i);
    if c ne 0 then printf "  [x^%o] %o\n", i, c; end if;
end for;
printf "\n";

// -----------------------------------------------------------------------
// Reduction: iterated Mobius search over a wide generator set.
// -----------------------------------------------------------------------
// Wide generator set for O_K: a + b*s5, a + b*w (w=(1+s5)/2), and small rationals.
w5 := (1+s5)/2;
B := 4;  // box bound
gens := [];
for a in [-B..B] do for b in [-B..B] do
    Append(~gens, k!a + k!b*k.1);
    Append(~gens, k!a + k!b*w5);
    Append(~gens, k!a/2 + k!b*k.1);
    Append(~gens, k!a + k!b*k.1/2);
end for; end for;
gens := Setseq(Seqset(gens));
printf "Generator set: %o elements\n\n", #gens;

best_f := f5; best_h := LogHeight(f5, k);
n5 := Degree(f5);

for pass in [1..3] do
    printf "=== Mobius pass %o (current h=%.4o) ===\n", pass, best_h;
    improved := false;
    for a in gens do for b in gens do for c in gens do for d in gens do
        det := a*d - b*c;
        if det eq 0 then continue; end if;
        g := MobiusPoly(best_f, a, b, c, d, Rk);
        if Degree(g) ne n5 then continue; end if;
        if LeadingCoefficient(g) eq 0 then continue; end if;
        h := LogHeight(g, k);
        if h lt best_h - 0.05 then
            best_h := h; best_f := g; improved := true;
            printf "  pass %o: (%o,%o,%o,%o) h=%.4o\n", pass, a,b,c,d, h;
        end if;
    end for; end for; end for; end for;
    if not improved then
        printf "  No further improvement.\n\n";
        break;
    end if;
    printf "\n";
end for;

// -----------------------------------------------------------------------
// Print the best degree-5 model found.
// -----------------------------------------------------------------------
printf "\n=== Best degree-5 model found: log-height %.4o ===\n", best_h;
printf "Coefficients of y^2 = f(x):\n";
for i in [Degree(best_f)..0 by -1] do
    c := Coefficient(best_f, i);
    if c ne 0 then printf "  [x^%o] %o\n", i, c; end if;
end for;

// Also try: does this model have integral coefficients? Clear denominators.
printf "\nChecking for integral model...\n";
f_best := best_f;
lcm_den := 1;
for c in Coefficients(f_best) do
    for e in Eltseq(k!c) do
        lcm_den := LCM(lcm_den, Denominator(e));
    end for;
end for;
if lcm_den gt 1 then
    printf "  Denominators present (LCM=%o). Scaling y -> y/%o to clear:\n", lcm_den, lcm_den;
    f_int := Rk![k!c * lcm_den^2 : c in Coefficients(f_best)];
    printf "  After scaling: log-height %.4o\n", LogHeight(f_int, k);
    printf "  Coefficients:\n";
    for i in [Degree(f_int)..0 by -1] do
        c := Coefficient(f_int, i);
        if c ne 0 then printf "    [x^%o] %o\n", i, c; end if;
    end for;
else
    printf "  Model already has integral coefficients.\n";
end if;
