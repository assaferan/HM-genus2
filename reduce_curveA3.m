// reduce_curveA3.m -- Stoll-style optimal affine reduction for Curve A.
//
// Strategy:
//   1. Start from best-so-far polynomial (Mobius pass from reduce_curveA.m).
//   2. Compute the 6 roots under each real embedding of K=Q(sqrt(5)).
//   3. For each embedding, the optimal affine x -> alpha*x + beta is:
//        beta  = centroid of roots (shifts roots to mean 0)
//        alpha = 1/RMS(|roots - centroid|)  (scales root spread to ~1)
//   4. Combine the two embedding-level optimal parameters to a single
//      (alpha,beta) in K by algebraizing the numerical values.
//   5. Apply, iterate until the height stabilises.
//
// Run: magma reduce_curveA3.m

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

// Embed polynomial f (over k) under a real embedding s5 -> ev.
function EmbedPoly(f, k, ev, RR)
    Rrx<t> := PolynomialRing(RR);
    return Rrx![RR!Eltseq(k!c)[1] + RR!Eltseq(k!c)[2]*ev : c in Coefficients(f)];
end function;

// Given K-element q ≈ val (numerically), try to find a good approximation
// as a + b * s5 with a, b ∈ {-B/N .. B/N} for small N, B.
function BestOKApprox(val, sqr5, k, B, Dens)
    best := k!0; bestd := RealField(50)!10^50;
    for N in Dens do
        for a in [-B*N..B*N] do
            for b in [-B*N..B*N] do
                cand := k!(a+0*s5)/N + k!(0+b*s5)/N;  // need exact k-coercion
                // Re-evaluate as (a/N) + (b/N)*s5
                v := RealField(50)!a/N + (RealField(50)!b/N)*sqr5;
                dv := Abs(val - v);
                if dv lt bestd then bestd := dv; best := k!a/N + k!b/N*k.1; end if;
            end for;
        end for;
    end for;
    return best, bestd;
end function;

// -----------------------------------------------------------------------
// Build starting polynomial (from reduce_curveA.m best result)
// -----------------------------------------------------------------------
printf "=== Setup: reconstruct Cmin and apply best first-pass Mobius ===\n";
C := Genus2CurveFromIgusa(QI, k);
Cmin, d := MinimalTwist(C, k);
f0, _ := HyperellipticPolynomials(Cmin);
f_cur := MobiusPoly(f0, 1, -s5+2, 0, s5+2, Rk);
printf "Starting polynomial log-height: %.4o  (from first-pass Mobius)\n\n", LogHeight(f_cur, k);

// -----------------------------------------------------------------------
// Stoll-style iteration: find optimal affine transform at each step.
// -----------------------------------------------------------------------
sqr5_val := [Sqrt(RealField(50)!5), -Sqrt(RealField(50)!5)];
CC := ComplexField(100);

function OptimalAffineForEmb(f_emb, RR)
    // Given embedded polynomial, return (alpha, beta) such that
    // roots are centred and spread ~1.
    CC := ComplexField(Precision(RR) + 20);
    rts := [r[1] : r in Roots(f_emb, CC)];
    n := #rts;
    if n eq 0 then return RR!1, RR!0; end if;
    mu := &+[Real(r) : r in rts] / n;    // centroid (real part only — totally real curve)
    rts2 := [r - mu : r in rts];
    sigma := Sqrt(&+[Abs(r)^2 : r in rts2] / n);
    if sigma lt RR!0.0001 then sigma := RR!1; end if;
    return 1/sigma, -mu/sigma;  // alpha, beta such that x -> alpha*(x+beta) normalizes
end function;

best_f := f_cur; best_h := LogHeight(f_cur, k);
for iter in [1..5] do
    printf "=== Iteration %o (current h=%.4o) ===\n", iter, best_h;

    // Compute optimal (alpha,beta) for each embedding.
    alphas := []; betas := [];
    for ev in sqr5_val do
        f_emb := EmbedPoly(best_f, k, ev, RealField(50));
        al, be := OptimalAffineForEmb(f_emb, RealField(50));
        Append(~alphas, al); Append(~betas, be);
        printf "  emb s5->%.3o: optimal alpha=%.4o beta=%.4o\n",
            RealField(5)!ev, al, be;
    end for;

    // Find best O_K approximations to alpha and beta.
    // Use the mean of the two embeddings as target for a first guess.
    al_mean := (alphas[1] + alphas[2]) / 2;
    be_mean := (betas[1] + betas[2]) / 2;
    printf "  Mean alpha=%.4o  Mean beta=%.4o\n", al_mean, be_mean;

    // Search for a + b*s5 with small |a|,|b| fitting alphas[1] under embedding 1.
    // We try many candidates and pick the one that gives the best resulting height.
    dens := [1,2,3,4,5,6,7,8,9,10];
    BBox := 20;
    best_sub_h := best_h; best_alpha := k!1; best_beta := k!0;

    // Candidates for alpha: approximate alphas[1] (dominant embedding) as (a+b*s5)/N.
    printf "  Searching for best affine transform...\n";
    alpha_cands := [];
    for N in dens do
        a0 := Round(N * RealField(5)!al_mean);
        for da in [-3..3] do
            aa := a0 + da;
            // Estimate b from embedding 1: aa/N + b/N * sqr5[1] ≈ alphas[1]
            b0 := Round(N * (RealField(5)!alphas[1] - RealField(5)!aa/N) / RealField(5)!sqr5_val[1]);
            for db in [-3..3] do
                bb := b0 + db;
                if aa eq 0 and bb eq 0 then continue; end if;
                cand := k!aa/N + k!bb/N*k.1;
                if cand ne 0 then Append(~alpha_cands, cand); end if;
            end for;
        end for;
    end for;
    // Also try pure-rational alpha candidates
    for N in dens do
        for a in [-BBox..BBox] do
            cand := k!a/N;
            if cand ne 0 then Append(~alpha_cands, cand); end if;
        end for;
    end for;
    alpha_cands := Setseq(Seqset(alpha_cands));   // deduplicate

    beta_cands := [];
    for N in dens do
        b0 := Round(N * RealField(5)!be_mean);
        for db in [-5..5] do
            bb := b0 + db;
            // b0+bb rational + k.1 component from embedding estimate
            brat_approx := RealField(5)!be_mean - RealField(5)!bb/N;
            bk := Round(N * brat_approx / RealField(5)!sqr5_val[1]);
            for dbk in [-3..3] do
                cand := k!bb/N + k!(bk+dbk)/N*k.1;
                Append(~beta_cands, cand);
            end for;
        end for;
    end for;
    // Also try small integer/unit shifts
    units := [k| 0, 1, -1, s5, -s5, (1+s5)/2, (1-s5)/2, (-1+s5)/2, (-1-s5)/2];
    for u in units do Append(~beta_cands, u); end for;
    beta_cands := Setseq(Seqset(beta_cands));

    printf "  alpha_cands: %o,  beta_cands: %o\n", #alpha_cands, #beta_cands;

    for alpha in alpha_cands do
        for beta in beta_cands do
            // Transform x -> alpha*(x + beta) = alpha*x + alpha*beta
            // As Mobius (a=alpha, b=alpha*beta, c=0, d=1):
            g := MobiusPoly(best_f, alpha, alpha*beta, 0, k!1, Rk);
            if Degree(g) lt Degree(best_f) then continue; end if;
            if LeadingCoefficient(g) eq 0 then continue; end if;
            h := LogHeight(g, k);
            if h lt best_sub_h then
                best_sub_h := h; best_alpha := alpha; best_beta := beta;
                printf "  alpha=%o  beta=%o  h=%.4o\n", alpha, beta, h;
            end if;
        end for;
    end for;

    if best_sub_h lt best_h - 0.01 then
        best_f := MobiusPoly(best_f, best_alpha, best_alpha*best_beta, 0, k!1, Rk);
        best_h := best_sub_h;
        printf "  Improved: h=%.4o  (alpha=%o, beta=%o)\n\n", best_h, best_alpha, best_beta;
    else
        printf "  No further improvement in this iteration.\n\n";
        break;
    end if;
end for;

// -----------------------------------------------------------------------
// Also compare with the degree-5 model (apply same reduction to it).
// -----------------------------------------------------------------------
printf "=== Degree-5 model reduction ===\n";
C := Genus2CurveFromIgusa(QI, k);
Cmin, _ := MinimalTwist(C, k);
has5, C5 := HasOddDegreeModel(Cmin);
if has5 then
    f5, _ := HyperellipticPolynomials(C5);
    printf "Degree-5 model log-height: %.4o\n", LogHeight(f5, k);

    // Apply the same Stoll-style reduction to f5
    f5_cur := f5; best_h5 := LogHeight(f5, k);
    for iter in [1..3] do
        best_sub_h5 := best_h5;
        for ev_idx in [1..2] do
            f_emb := EmbedPoly(f5_cur, k, sqr5_val[ev_idx], RealField(50));
            al5, be5 := OptimalAffineForEmb(f_emb, RealField(50));
            // Try rational approximations of al5, be5
            for Na in [1..10] do
                aa := Round(Na * RealField(5)!al5);
                for da in [-2..2] do
                    alpha_try := k!(aa+da)/Na;
                    if alpha_try eq 0 then continue; end if;
                    mu := RealField(5)!be5;
                    for mb in [-10..10] do
                        beta_try := k!mb/Na;
                        g := MobiusPoly(f5_cur, alpha_try, alpha_try*beta_try, 0, k!1, Rk);
                        if Degree(g) lt Degree(f5_cur) or LeadingCoefficient(g) eq 0 then continue; end if;
                        h := LogHeight(g, k);
                        if h lt best_sub_h5 then
                            best_sub_h5 := h;
                            f5_cur := g;
                            printf "  deg5 iter %o: alpha=%o beta=%o h=%.4o\n", iter, alpha_try, beta_try, h;
                        end if;
                    end for;
                end for;
            end for;
        end for;
        if best_sub_h5 ge best_h5 - 0.01 then break; end if;
        best_h5 := best_sub_h5;
    end for;
    printf "Best degree-5 model h=%.4o\n", best_h5;
else
    printf "No degree-5 model.\n";
end if;

// -----------------------------------------------------------------------
// Summary
// -----------------------------------------------------------------------
printf "\n=== Final summary ===\n";
printf "Degree-6 model log-height: %.4o\n", best_h;
printf "y^2 = f(x), f =\n";
for i in [Degree(best_f)..0 by -1] do
    c := Coefficient(best_f, i);
    if c ne 0 then printf "  [x^%o] %o\n", i, c; end if;
end for;
