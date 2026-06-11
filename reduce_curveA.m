// reduce_curveA.m -- Search for small-coefficient models of Curve A over K=Q(sqrt(5)).
//
// Approach:
//   0. Print the raw Cmin polynomial from MinimalTwist.
//   1. O_K shifts: x -> x + alpha, alpha = a + b*s5, |a|,|b| <= ShiftBox.
//   2. Flip + shift: x -> 1/x then shift.
//   3. Full Mobius search: x -> (a*x+b)/(c*x+d) with generators from a small O_K set.
//   4. Stoll-style reduction (one step): embed under both real places, find the
//      GL2(R) transform that equidistributes the roots, approximate by GL2(O_K).
//      This is the most principled approach and can be iterated.
//
// Run: magma reduce_curveA.m

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

// -----------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------

// Log-height of a polynomial over k: max of log|coefficients| under both embeddings.
function LogHeight(f, k)
    sqr5 := Sqrt(RealField(30)!5);
    embs := [sqr5, -sqr5];
    h := RealField(30)!0;
    for c in Coefficients(f) do
        seq := Eltseq(k!c);
        a := RealField(30)!seq[1]; b := RealField(30)!seq[2];
        for e in embs do
            h := Max(h, Abs(a + b*e));
        end for;
    end for;
    return h gt 0 select Log(h) else RealField(30)!0;
end function;

// Apply Mobius x -> (a*x+b)/(c*x+d) to a degree-n poly (returns numerator poly).
function MobiusPoly(f, a, b, c, d, Rk)
    n := Degree(f); g := Rk!0;
    for i in [0..n] do
        g +:= Coefficient(f,i) * (a*Rk.1 + b)^i * (c*Rk.1 + d)^(n-i);
    end for;
    return g;
end function;

// Remove leading rational scalar (divide by leading coefficient over k).
function Normalise(f, k)
    lc := LeadingCoefficient(f);
    if lc eq 0 then return f; end if;
    return f / lc;
end function;

// Print the polynomial neatly.
procedure PrintPoly(f, k)
    n := Degree(f);
    for i in [n..0 by -1] do
        c := Coefficient(f, i);
        if c ne 0 then printf "  [x^%o] %o\n", i, c; end if;
    end for;
end procedure;

// -----------------------------------------------------------------------
// Step 0: Reconstruct Cmin and print its polynomial.
// -----------------------------------------------------------------------
printf "=== Step 0: Cmin polynomial ===\n";
C := Genus2CurveFromIgusa(QI, k);
Cmin, d := MinimalTwist(C, k);
f0, _ := HyperellipticPolynomials(Cmin);
printf "Twist d = %o\n", d;
printf "f0 (degree %o, log-height %.2o):\n", Degree(f0), LogHeight(f0, k);
PrintPoly(f0, k);
printf "\n";

best_h := LogHeight(f0, k); best_f := f0; best_desc := "Cmin";

// -----------------------------------------------------------------------
// Step 1: O_K shifts   x -> x + alpha,  |a|, |b| <= ShiftBox
// -----------------------------------------------------------------------
ShiftBox := 10;
printf "=== Step 1: O_K shifts (a+b*s5, |a|,|b| <= %o) ===\n", ShiftBox;
for a in [-ShiftBox..ShiftBox] do
    for b in [-ShiftBox..ShiftBox] do
        alpha := k!a + k!b * k.1;
        fa := Evaluate(f0, x + alpha);
        h := LogHeight(fa, k);
        if h lt best_h then
            best_h := h; best_f := fa;
            best_desc := Sprintf("shift x+(%o+%o*s5)", a, b);
            printf "  a=%o+%o*s5  h=%.4o\n", a, b, h;
        end if;
    end for;
end for;
printf "  Best so far: %o  h=%.4o\n\n", best_desc, best_h;

// -----------------------------------------------------------------------
// Step 2: Flip + O_K shift
// -----------------------------------------------------------------------
printf "=== Step 2: Flip + O_K shifts ===\n";
f_flip := Polynomial([Coefficient(f0, Degree(f0)-i) : i in [0..Degree(f0)]]);
for a in [-ShiftBox..ShiftBox] do
    for b in [-ShiftBox..ShiftBox] do
        alpha := k!a + k!b * k.1;
        fa := Evaluate(f_flip, x + alpha);
        h := LogHeight(fa, k);
        if h lt best_h then
            best_h := h; best_f := fa;
            best_desc := Sprintf("flip+shift(%o+%o*s5)", a, b);
            printf "  flip+a=%o+%o*s5  h=%.4o\n", a, b, h;
        end if;
    end for;
end for;
printf "  Best so far: %o  h=%.4o\n\n", best_desc, best_h;

// -----------------------------------------------------------------------
// Step 3: Mobius search  x -> (a*x+b)/(c*x+d),  a,b,c,d from small O_K set
// -----------------------------------------------------------------------
printf "=== Step 3: Mobius search ===\n";
// Generators to try for each entry (elements of O_K with small norm)
gens := [k| 0, 1, -1, s5, -s5, 1+s5, 1-s5, -1+s5, -1-s5, 2, -2, 2+s5, 2-s5, -2+s5, -2-s5];
n6 := Degree(f0);
nmob := 0; ncheck := 0;
for a in gens do for b in gens do for c in gens do for d in gens do
    det := a*d - b*c;
    if det eq 0 then continue; end if;
    nmob +:= 1;
    g := MobiusPoly(f0, a, b, c, d, Rk);
    if Degree(g) lt n6 then continue; end if;
    lc := LeadingCoefficient(g);
    if lc eq 0 then continue; end if;
    ncheck +:= 1;
    h := LogHeight(g, k);
    if h lt best_h then
        best_h := h; best_f := g;
        best_desc := Sprintf("Mobius (%o,%o,%o,%o)", a, b, c, d);
        printf "  (%o,%o,%o,%o)  h=%.4o\n", a,b,c,d,h;
    end if;
end for; end for; end for; end for;
printf "  Tried %o transforms (%o non-degenerate). Best: %o  h=%.4o\n\n",
    nmob, ncheck, best_desc, best_h;

// -----------------------------------------------------------------------
// Step 4: Stoll-style one-step reduction via real embeddings.
//
// For K = Q(sqrt(5)) (totally real), embed f under both places:
//   sigma_1: s5 -> +sqrt(5)   sigma_2: s5 -> -sqrt(5)
// Find the complex roots of each embedded poly, compute the GL2(R) Minkowski
// reduction matrix for each embedding, average them, round to GL2(O_K),
// and apply the resulting Mobius transform.
// -----------------------------------------------------------------------
printf "=== Step 4: Stoll-style reduction via real embeddings ===\n";

function EmbedPoly(f, k, sqr5_val, RR)
    // Embed f over k -> polynomial over RR via s5 -> sqr5_val.
    Rrx<t> := PolynomialRing(RR);
    return Rrx![RR!Eltseq(k!c)[1] + RR!Eltseq(k!c)[2]*sqr5_val : c in Coefficients(f)];
end function;

function StollMatrixFromRoots(f_real, RR)
    // Given a real polynomial, compute GL2(R) that "centres" the roots.
    // Strategy: find the centroid of roots (as Mobius input), translate to 0,
    // then scale to unit diameter.
    CC := ComplexField(Precision(RR));
    rts := [r[1] : r in Roots(f_real, CC)];
    if #rts eq 0 then return IdentityMatrix(RR, 2); end if;
    mu := &+rts / #rts;     // centroid
    re_mu := Real(mu);
    rts2 := [r - mu : r in rts];
    sigma := Sqrt(&+[Abs(r)^2 : r in rts2] / #rts2 + RR!1);
    // M = [[1/sigma, -re_mu/sigma],[0, 1]] (translate then scale)
    return Matrix(RR, 2, 2, [1/sigma, -re_mu/sigma, 0, 1]);
end function;

RR := RealField(50);
sqr5 := Sqrt(RR!5);
emb_vals := [sqr5, -sqr5];
Ms := [];
for ev in emb_vals do
    f_emb := EmbedPoly(f0, k, ev, RR);
    M := StollMatrixFromRoots(f_emb, RR);
    Append(~Ms, M);
end for;

// Average the two matrices, then round to GL2(Z) (since for small search
// the off-diagonal might be Z-approximable; if it has a sqrt(5) component
// we'll catch it in the Mobius search above).
M_avg := (Ms[1] + Ms[2]) / 2;
printf "  Average GL2(R) matrix:\n    [%.4o  %.4o]\n    [%.4o  %.4o]\n",
    M_avg[1][1], M_avg[1][2], M_avg[2][1], M_avg[2][2];

// Try rounding to O_K entries: approximate each entry as a+b*s5 with small |a|,|b|.
function RoundToOK(r, k, sqr5, B)
    // Approximate r ∈ R as a+b*s5 with a,b ∈ {-B..B}/N for denominators N=1,2,3,4.
    best := [Rationals()|0, 0]; bestd := Abs(r);
    for N in [1,2,3,4] do
        for a in [-B*N..B*N] do
            for b in [-B*N..B*N] do
                v := a/N + (b/N)*sqr5;
                dv := Abs(r - v);
                if dv lt bestd then bestd := dv; best := [a/N, b/N]; end if;
            end for;
        end for;
    end for;
    return k!best[1] + k!best[2]*k.1;
end function;

// Round entries of M_avg to O_K (or small rational multiples).
a_r := RoundToOK(M_avg[1][1], k, sqr5, 3);
b_r := RoundToOK(M_avg[1][2], k, sqr5, 3);
c_r := RoundToOK(M_avg[2][1], k, sqr5, 3);
d_r := RoundToOK(M_avg[2][2], k, sqr5, 3);
printf "  Rounded to O_K: a=%o b=%o c=%o d=%o\n", a_r, b_r, c_r, d_r;

det_r := a_r*d_r - b_r*c_r;
if det_r ne 0 then
    g_stoll := MobiusPoly(f0, a_r, b_r, c_r, d_r, Rk);
    if Degree(g_stoll) eq n6 then
        h_stoll := LogHeight(g_stoll, k);
        printf "  Stoll-step transform: h=%.4o\n", h_stoll;
        if h_stoll lt best_h then
            best_h := h_stoll; best_f := g_stoll;
            best_desc := Sprintf("Stoll-step (%o,%o,%o,%o)", a_r,b_r,c_r,d_r);
        end if;
    else
        printf "  Stoll-step: degree dropped (degenerate).\n";
    end if;
else
    printf "  Stoll-step: rounded matrix singular.\n";
end if;
printf "\n";

// -----------------------------------------------------------------------
// Summary
// -----------------------------------------------------------------------
printf "=== Best model: %o  (log-height %.4o) ===\n", best_desc, best_h;
printf "y^2 = f(x), f =\n";
PrintPoly(best_f, k);
