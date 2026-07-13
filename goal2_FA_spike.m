// goal2_FA_spike.m -- feasibility spike for the degree-400 F_A resolvent (idx 33).
//
// F_A(proj, F_7-lines): the field of the 400 lines of P^3(F_7) in A[7].  It needs no
// arithmetic F_49-structure J (unlike the clean degree-50 field), so it IS formable
// from the numerical 7-torsion.  This spike MEASURES whether it is recognizable:
//   (1) time per torsion point (-> x2400 total),
//   (2) the coefficient height of R(X) = prod_lines (X - theta_line), which sets the
//       precision floor for recognizing R(X) in Q[X].
//
// Kummer coordinate: FromAnalyticJacobian gives the degree-2 divisor {P1,P2}; its
// Mumford u-poly u(x)=x^2 - s1 x + s2 (s1=x1+x2, s2=x1 x2) is +-1-invariant.
// Line invariant: theta_line = sum over the 3 Kummer points {[P],[2P],[3P]} of
// t = s1 + rho*s2  (rho a fixed rational), which is invariant under the F_7^* action
// and Galois-equivariant, so the multiset {theta_line} is Gal(Qbar/Q)-stable and
// R(X) in Q[X], degree 400.
//
// Run: magma goal2_FA_spike.m
SetColumns(0);
prec := 80;
rho  := 3;                       // fixed rational mixing s1,s2 (separates lines)
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
C := HyperellipticCurve(-10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9);
LOGF := "goal2_FA_spike_out.txt";
PrintFile(LOGF, "# F_A degree-400 resolvent feasibility spike (idx 33)" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;
LOG(Sprintf("precision = %o digits, rho = %o", prec, rho));

t0 := Cputime();
AJ  := AnalyticJacobian(C : Precision := prec);
BPM := BigPeriodMatrix(AJ);          // 2 x 4
CC  := BaseRing(BPM);
lat := [ Matrix(CC,2,1,[BPM[1,j],BPM[2,j]]) : j in [1..4] ];
LOG(Sprintf("period matrix built (%.1o s)", Cputime(t0)));

// Kummer coordinate t of the torsion point with F_7-vector w = (w1..w4) in {0..6}^4.
// returns <ok, t> ; ok=false if the divisor is anomalous (point at infinity etc.)
function KummerT(w)
    z := (CC!0) * lat[1];
    for j in [1..4] do z +:= (CC!w[j]) * lat[j]; end for;
    z := (1/7) * z;
    pts := FromAnalyticJacobian(z, AJ);
    if #pts ne 2 then return false, CC!0; end if;   // anomalous divisor (pt at infinity)
    x1 := pts[1][1]; x2 := pts[2][1];               // affine x-coordinates <x,y>
    s1 := x1 + x2;                                   // Mumford u(x) = x^2 - s1 x + s2
    s2 := x1 * x2;
    return true, s1 + rho*s2;
end function;

// ---- enumerate the 400 lines of P^3(F_7): vectors normalized to leading entry 1 ----
lines := [];
for v1 in [0..6] do for v2 in [0..6] do for v3 in [0..6] do for v4 in [0..6] do
    w := [v1,v2,v3,v4];
    // leading (first) nonzero entry must be 1  ->  one representative per line
    lead := 0;
    for j in [1..4] do if w[j] ne 0 then lead := w[j]; break; end if; end for;
    if lead eq 1 then Append(~lines, w); end if;
end for; end for; end for; end for;
LOG(Sprintf("enumerated %o lines of P^3(F_7)", #lines));

// ---- timing probe on the first few Kummer evaluations ----
tp := Cputime();
nprobe := 12; got := 0;
for i in [1..nprobe] do
    ok, _ := KummerT(lines[i]);
    if ok then got +:= 1; end if;
end for;
per := Cputime(tp)/nprobe;
LOG(Sprintf("per-Kummer-eval time ~ %.3o s  (%o/%o clean on probe)", per, got, nprobe));
LOG(Sprintf("=> full 400 lines x 3 Kummer pts = 1200 evals ~ %.1o s projected", per*1200));

// ---- compute theta_line for all 400 lines ----
tt := Cputime();
thetas := []; clean := 0; skipped := 0;
for li in [1..#lines] do
    w := lines[li];
    // the 3 Kummer points of the line: k*w mod 7, k=1,2,3
    vals := []; okline := true;
    for k in [1..3] do
        wk := [ (k*w[j]) mod 7 : j in [1..4] ];
        ok, t := KummerT(wk);
        if not ok then okline := false; break; end if;
        Append(~vals, t);
    end for;
    if okline then
        Append(~thetas, vals[1]+vals[2]+vals[3]);
        clean +:= 1;
    else
        skipped +:= 1;
    end if;
    if li mod 50 eq 0 then
        LOG(Sprintf("  ... %o/%o lines (%o clean, %o skipped, %.0o s)", li, #lines, clean, skipped, Cputime(tt)));
    end if;
end for;
LOG(Sprintf("theta_line computed: %o clean, %o skipped (%.1o s)", clean, skipped, Cputime(tt)));

// typical magnitude of theta_line
mags := [ Log(10, Abs(t)+10^(-prec)) : t in thetas | Abs(t) gt 0 ];
avgmag := &+mags / #mags;
maxmag := Max(mags);
LOG(Sprintf("theta_line magnitude: mean 10^%.2o, max 10^%.2o (n=%o)", avgmag, maxmag, #thetas));

// ---- build the resolvent (numerically) over the clean lines and read coefficient heights ----
CCx<X> := PolynomialRing(CC);
R := &*[ X - t : t in thetas ];
coeffs := Coefficients(R);
cmags := [ Log(10, Abs(c)+10^(-prec)) : c in coeffs | Abs(c) gt 10^(-prec) ];
LOG(Sprintf("\nresolvent degree = %o (clean lines)", #thetas));
LOG(Sprintf("max |coefficient| ~ 10^%.1o  <-- precision floor to even REPRESENT it", Max(cmags)));
LOG(Sprintf("(theory: a degree-400 product peaks near 10^(120 + 200*meanmag) = 10^%.0o)", 120 + 200*avgmag));

// try to recognize the trace and product (cheapest coefficients) as rationals
tr := -Coefficient(R, #thetas-1);         // sum of theta_line
LOG(Sprintf("\ntrace (sum theta_line) ~ %.6o + %.6o i", Real(tr), Imaginary(tr)));
LOG(Sprintf("  |Im(trace)| = 10^%.1o (should be ~0 if genuinely rational)", Log(10, Abs(Imaginary(tr))+10^(-prec))));
try
    q := BestApproximation(Real(tr), 10^(prec div 2));
    LOG(Sprintf("  trace best rational approx: %o", q));
catch e; LOG("  trace recognition failed"); end try;

LOG("\n=== VERDICT ===");
LOG(Sprintf("degree 400, max coeff ~ 10^%.0o.  Recognizing R(X) in Q[X] needs precision", Max(cmags)));
LOG(Sprintf(">~ %o digits (plus denominator height at {2,3,5,7}).  Compare: prec %o here.", Ceiling(Max(cmags)), prec));
LOG("DONE.");
exit;
