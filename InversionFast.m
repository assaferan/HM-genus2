// Optimized inverse map x -> tau (production version; see Inversion.m for the
// reference). Three speedups over the reference:
//   1. ThetaJet computes the theta values AND their tau-derivatives in a single
//      lattice pass -> an ANALYTIC Jacobian (no finite differencing).
//   2. FindTauFast uses Levenberg-Marquardt (damped Gauss-Newton) with
//      backtracking and positive-definiteness guards.
//   3. The lattice truncation N is adaptive: coarse while the residual is large,
//      ramped to full precision only near the solution.
import "HMSurface.m": HMSurface;
import "Inversion.m": ThetaPt, QEval, Residual, verify_tau;

// prec := 200;                       // push it freely now
// CC := ComplexField(prec);
//PI := Pi(CC); ii := CC.1;
//P5<x0,x1,x2,x3,x4> := PolynomialRing(CC, 5);
// qs := [...];  // quintics (+ a few sextics), coefficients coerced K -> CC
// dqs := [ [ Derivative(q, l) : l in [1..5] ] : q in qs ];   // partials, once

// Theta values + their three tau-derivatives (d/da, d/db, d/dc) for all five
// theta coordinates, computed in ONE lattice pass. Each term's tau-derivative
// just multiplies the term E by pi*i times the relevant quadratic monomial in v,
// so the derivative sums piggyback on the value sum at negligible extra cost.
function ThetaJet(z, a, b, c, N, CC)
    PI := Pi(CC); ii := CC.1;
    th := [CC|0,0,0,0,0]; tha := th; thb := th; thc := th;
    for j in [0..4] do
        sj:=CC!0; sa:=CC!0; sb:=CC!0; sc:=CC!0;
        for n1 in [-N..N], n2 in [-N..N] do
            v1 := CC!n1; v2 := n2 + CC!j/5;
            E := Exp(PI*ii*(a*v1^2 + 2*b*v1*v2 + c*v2^2)
                     + 2*PI*ii*(v1*z[1] + v2*z[2]));
            sj +:= E; sa +:= (PI*ii*v1^2)*E;            // d/da: factor pi i v1^2
            sb +:= (PI*ii*2*v1*v2)*E; sc +:= (PI*ii*v2^2)*E;  // d/db, d/dc
        end for;
        th[j+1]:=sj; tha[j+1]:=sa; thb[j+1]:=sb; thc[j+1]:=sc;
    end for;
    return th, tha, thb, thc;
end function;

// Residual vector R and its analytic Jacobian J (w.r.t. a,b,c) in one shot.
// dqs are the precomputed partials of the quintics. For each z we dehomogenize
// the theta point in its largest coordinate (a numerically safe affine chart),
// push the tau-derivatives through the quotient rule, and chain-rule with the
// quintic gradients to get d/d{a,b,c} of each residual entry.
function ResJac(qs, dqs, zs, a, b, c, N, CC)
    Rseq:=[CC|]; JA:=[CC|]; JB:=[CC|]; JC:=[CC|];
    for z in zs do
        th,tha,thb,thc := ThetaJet(z,a,b,c,N, CC);
        _, k := Max([Abs(th[l]) : l in [1..5]]);     // largest coord = safe chart
        mu:=th[k]; ma:=tha[k]; mb:=thb[k]; mc:=thc[k];
        Pt := [ th[l]/mu : l in [1..5] ];                     // dehomogenized point
        Pa := [ (tha[l]*mu - th[l]*ma)/mu^2 : l in [1..5] ];  // d Pt / da (quotient rule)
        Pb := [ (thb[l]*mu - th[l]*mb)/mu^2 : l in [1..5] ];  // d Pt / db
        Pc := [ (thc[l]*mu - th[l]*mc)/mu^2 : l in [1..5] ];  // d Pt / dc
        for i in [1..#qs] do
            Append(~Rseq, Evaluate(qs[i], Pt));
            grad := [ Evaluate(dqs[i][l], Pt) : l in [1..5] ];   // grad q_i at Pt
            // chain rule: d q_i / d{a,b,c} = grad q_i . d Pt / d{a,b,c}
            Append(~JA, &+[ grad[l]*Pa[l] : l in [1..5] ]);
            Append(~JB, &+[ grad[l]*Pb[l] : l in [1..5] ]);
            Append(~JC, &+[ grad[l]*Pc[l] : l in [1..5] ]);
        end for;
    end for;
    m := #Rseq;
    R := Matrix(CC, m, 1, Rseq);
    J := ZeroMatrix(CC, m, 3);
    for r in [1..m] do J[r,1]:=JA[r]; J[r,2]:=JB[r]; J[r,3]:=JC[r]; end for;
    return R, J;
end function;

/* 
function FindTau(qs, dqs, tau0, zs, N : iters := 40)
    a:=tau0[1]; b:=tau0[2]; c:=tau0[3];
    for it in [1..iters] do
        R, J := ResJac(qs,dqs,zs,a,b,c,N); m := Nrows(R);
        Js := ZeroMatrix(CC,3,m);
        for r in [1..m], kk in [1..3] do Js[kk,r] := Conjugate(J[r,kk]); end for;
        if Rank(Js*J) lt 3 then break; end if;
        del := -(Js*J)^(-1) * (Js*R);
        a +:= del[1,1]; b +:= del[2,1]; c +:= del[3,1];
        nrm := Sqrt(&+[Abs(R[r,1])^2 : r in [1..m]]);
        printf "it %o: |R| = %o\n", it, RealField(8)!nrm;
        if nrm lt 10^(-(prec-10)) then break; end if;
    end for;
    return [a,b,c];
end function; 
*/

// Levenberg-Marquardt solve for tau, starting from tau0. lambda is the damping:
// small -> Gauss-Newton (fast near solution), large -> gradient descent (safe far
// away). Each accepted step must keep Im(tau) positive definite and strictly
// decrease the residual norm (backtracking by inflating lambda otherwise).
// stag_floor: only abandon while |R| is still this large (we are far from a zero).
// stag_rel/stag_win: if the relative improvement in |R| stays below stag_rel for
//   stag_win consecutive iterations, the trial has settled into a spurious fixed
//   point of the least-squares objective (not a zero) -- abandon it early instead
//   of burning the remaining iterations.
// lm_floor: abandon a trial whose Im(tau) smallest eigenvalue collapses below this.
//   Such trials are wandering into a degenerate region (theta convergence dies and
//   N explodes); the true basin keeps lm well away from 0, so this only prunes
//   expensive dead ends. Set to 0 to disable (e.g. for the in-basin polish).
// tol_pow: stop when |R| < 10^(-tol_pow). Default 0 => 3*prec/7 (loose; fine for
//   basin detection in the search). The polish needs tau accurate to nearly the
//   FULL working precision (ReconstructCurve's LLL fails on periods that are only
//   correct to a fraction of prec), so it passes tol_pow ~ prec - buffer.
function FindTauFast(qs, dqs, tau0, zs, CC : iters := 60, verbose := false,
                     stag_floor := 1e-3, stag_rel := 1e-3, stag_win := 4, coarse := true,
                     lm_floor := 1e-2, tol_pow := 0)
    prec := Precision(CC);
    a := tau0[1]; b := tau0[2]; c := tau0[3];
    lambda := CC!1e-4;                       // LM damping parameter
    nrm := 10^9;
    tol := tol_pow gt 0 select Real(CC!10^(-tol_pow))   // explicit tolerance
                          else Real(CC!10^(-Floor(3*prec/7)));  // default (loose)
    prevnrm := 10^18; stag := 0;             // stagnation tracking (prevnrm large to start)
    for it in [1..iters] do
        // Smallest eigenvalue of Im(tau): controls theta convergence rate.
        lm := ((Imaginary(a)+Imaginary(c))
               - Sqrt((Imaginary(a)-Imaginary(c))^2+4*Imaginary(b)^2))/2;
        if lm_floor gt 0 and lm lt lm_floor then
            if verbose then printf "  [Im(tau) eigenvalue %o below floor, abandoning trial]\n", RealField(6)!lm; end if;
            break;
        end if;
        // Truncation N needed so tail terms drop below the target precision.
        // lm is real (built from Im parts); keep the floor real so Max is ordered.
        N_need := Ceiling(Sqrt(prec*Log(10.0)/(Real(Pi(CC))*Max(lm, Real(CC!0.005)))))+3;
        // Coarse theta sums while the residual is large (each iteration is O(N^2)).
        // NB: with coarse N the residual cannot drop below the truncation floor,
        // so an in-basin polish (where we need full accuracy) must set coarse:=false
        // to use full N from the start -- otherwise it stalls at the floor.
        if not coarse then
            N := N_need;
        elif nrm gt 10^(-2) then
            N := Min(N_need, Max(8, N_need div 3));
        elif nrm gt 10^(-6) then
            N := Min(N_need, Max(12, (2*N_need) div 3));
        else
            N := N_need;
        end if;
        R, J := ResJac(qs, dqs, zs, a, b, c, N, CC);
        m := Nrows(R);
        nrm := Sqrt(&+[Abs(R[i,1])^2 : i in [1..m]]);
        if verbose then
            printf "it %o: N=%o lam_min=%o |R|=%o\n", it, N, RealField(6)!lm, RealField(8)!nrm;
        end if;
        if nrm lt tol then break; end if;

        // Stagnation check: abandon trials stuck at a spurious (nonzero) fixed point.
        if nrm gt stag_floor and nrm ge prevnrm*(1 - stag_rel) then
            stag +:= 1;
        else
            stag := 0;
        end if;
        if stag ge stag_win then
            if verbose then printf "  [stalled at |R|=%o, abandoning trial]\n", RealField(8)!nrm; end if;
            break;
        end if;
        prevnrm := nrm;

        Js := ZeroMatrix(CC, 3, m);              // conjugate transpose J*
        for i in [1..m], k in [1..3] do Js[k,i] := Conjugate(J[i,k]); end for;
        JsJ := Js*J;                             // 3x3 normal matrix (Hermitian, PSD)
        // Rank-deficient normal matrix: the residual fails to constrain all three
        // tau directions here (a degenerate start). Even heavy LM damping by D can
        // stay singular if a Jacobian column vanishes, so abandon this trial -- the
        // next random restart will sample elsewhere.
        if Rank(JsJ) lt 3 then
            if verbose then printf "  [singular normal matrix, abandoning trial]\n"; end if;
            break;
        end if;

        accepted := false;
        // Try the undamped Gauss-Newton step first (JsJ is now nonsingular).
        step := -(JsJ)^(-1)*(Js*R);
        an := a + step[1,1]; bn := b + step[2,1]; cn := c + step[3,1];
        pd := Imaginary(an) gt 0 and Imaginary(an)*Imaginary(cn) - Imaginary(bn)^2 gt 0;
        if pd then
            Rn, _ := ResJac(qs, dqs, zs, an, bn, cn, N, CC);
            nn := Sqrt(&+[Abs(Rn[i,1])^2 : i in [1..Nrows(Rn)]]);
            if nn lt nrm then
                a := an; b := bn; c := cn; lambda := CC!1e-4; accepted := true;
            end if;
        end if;

        // If GN was rejected, fall back to damped LM steps, inflating lambda
        // until a step is accepted (or lambda blows up -> give up this iter).
        D := DiagonalMatrix([JsJ[i,i] : i in [1..3]]);   // LM scaling (diag of J* J)
        while not accepted and Abs(lambda) lt 1e18 do
            step := -(JsJ + lambda*D)^(-1)*(Js*R);
            an := a + step[1,1]; bn := b + step[2,1]; cn := c + step[3,1];
            pd := Imaginary(an) gt 0 and Imaginary(an)*Imaginary(cn) - Imaginary(bn)^2 gt 0;
            if pd then
                Rn, _ := ResJac(qs, dqs, zs, an, bn, cn, N, CC);
                nn := Sqrt(&+[Abs(Rn[i,1])^2 : i in [1..Nrows(Rn)]]);
                if nn lt nrm then
                    a := an; b := bn; c := cn; lambda := lambda/4; accepted := true;
                else
                    lambda := lambda*3;
                end if;
            else
                lambda := lambda*3;
            end if;
        end while;
    end for;
    return [a,b,c], nrm;
end function;

// Complex embedding K -> CC for the base field K (Rationals or a number field).
// For Q it is the canonical inclusion; for a number field it is the Embedding-th
// root of the defining polynomial, with roots sorted canonically (by real then
// imaginary part) so the SAME embedding is chosen at every precision. Different
// embeddings give Galois-conjugate inputs (hence conjugate curves).
function ComplexEmbedding(K, idx, CC)
    if Type(K) eq FldRat then return hom< K -> CC | >; end if;
    rts := [ r[1] : r in Roots(ChangeRing(DefiningPolynomial(K), CC)) ];
    Sort(~rts, func< a, b | Sign(Real(a) - Real(b)) ne 0
                            select Sign(Real(a) - Real(b))
                            else   Sign(Imaginary(a) - Imaginary(b)) >);
    error if idx gt #rts, "ComplexEmbedding: Embedding index out of range; K has", #rts, "embeddings";
    return hom< K -> CC | rts[idx] >;
end function;

// Coerce the K-quintics qs into CC via the embedding emb; also return their partials.
function CoerceQuintics(qs, K, emb, CC)
    P5K := Parent(qs[1]);
    P5CC := PolynomialRing(CC, 5);
    if Type(K) eq FldRat then
        cm := hom< P5K -> P5CC | [P5CC.j : j in [1..5]] >;          // Q-coeffs coerce canonically
    else
        cm := hom< P5K -> P5CC | emb, [P5CC.j : j in [1..5]] >;      // via the chosen embedding
    end if;
    qsC := [ cm(q) : q in qs ];
    return qsC, [ [ Derivative(q, l) : l in [1..5] ] : q in qsC ];
end function;

// x in P^3(K) to tau such that A_tau = A_x (over the chosen embedding of K).
//
// Two-phase strategy (the residual surface has many spurious local minima, so a
// random start converges to the true basin only ~10% of the time):
//   (A) LOCATE -- cheap random restarts at reduced precision; stop at the FIRST
//       start that lands in the basin (score below basin_thresh). Doing the
//       restarts cheaply, and stopping early, is the dominant speedup.
//   (B) POLISH -- refine that single winner to full precision, then verify.
//
// K selects the base field (default Q). Over a number field only SOME complex
// embeddings yield a positive-definite period tau, and which ones is x-dependent;
// so Embedding := 0 (the default) tries all embeddings round-robin and keeps the
// first that lands in a basin, while Embedding := i forces the i-th embedding. The
// K-quintics are coerced to CC via the chosen embedding before the numeric search;
// the winning embedding is returned so the caller can recognize invariants in K.
//
// Only ~10% of random starts land in the true basin, so trials is large; the
// early break means we usually stop after ~10 (P(miss all 100) < 1e-4).
function x_to_tau_fast(x, CC : K := Rationals(), Embedding := 0, trials := 100,
                       search_prec := 0, verbose := false, basin_thresh := 1e-6, polish_iters := 40)
    prec := Precision(CC);
    // HM surface over K and its defining quintics (exact, precision-independent).
    P := ProjectiveSpace(K, 4);
    F := HorrocksMumfordBundle(P);
    A := HMSurface(F, x);
    qs := [b : b in MinimalBasis(Ideal(A)) | Degree(b) eq 5];   // the 3 defining quintics
    assert #qs eq 3;

    // ---- Phase A: locate the basin at reduced precision ----
    sp := search_prec gt 0 select search_prec else Max(30, Min(prec, 40));
    CCs := ComplexFieldExtra(sp); iis := CCs.1;
    nemb := Type(K) eq FldRat select 1 else Degree(K);
    embset := Embedding eq 0 select [1..nemb] else [Embedding];
    // Pre-coerce the quintics once per candidate embedding (at search precision).
    coercedS := [* *];
    for e in embset do
        qe, dqe := CoerceQuintics(qs, K, ComplexEmbedding(K, e, CCs), CCs);
        Append(~coercedS, <qe, dqe>);
    end for;
    zs := [ [CCs| (Random(-50,50)+iis*Random(10,40))/100,
                  (Random(-50,50)+iis*Random(10,40))/100 ] : kk in [1..6] ];
    Nscore := Max(14, Ceiling(Sqrt(sp/4)));       // truncation for scoring
    best := []; bestR := Infinity(); chosen := embset[1];
    for trial in [1..trials] do
        j := ((trial-1) mod #embset) + 1;          // round-robin over candidate embeddings
        qsS := coercedS[j][1]; dqsS := coercedS[j][2];
        // Wide sampling box: the true tau routinely lies well outside a reduced
        // box (e.g. for x=[1,2,3,4]: a=-1.26+0.24i, b=0.95-0.05i, c=-2.04+1.92i),
        // so a narrow box rarely seeds the basin. This roughly doubles the hit rate.
        tau0 := [ (Random(-200,200)+iis*Random(10,120))/100,   // a
                  (Random(-150,150)+iis*Random(-50,50))/100,   // b (now complex)
                  (Random(-250,250)+iis*Random(10,250))/100 ]; // c
        ts, _ := FindTauFast(qsS, dqsS, tau0, zs, CCs : iters := 50, verbose := verbose);
        // Accept only positive-definite Im(tau); keep the lowest-residual one.
        posdef := Imaginary(ts[1]) gt 0 and
                  Imaginary(ts[1])*Imaginary(ts[3]) - Imaginary(ts[2])^2 gt 0;
        r := Sqrt(&+[Abs(rv)^2 : rv in Residual(qsS, zs, ts[1], ts[2], ts[3], Nscore, CCs)]);
        if verbose then printf "search trial %o (embedding %o): posdef=%o |R|=%o\n", trial, embset[j], posdef, RealField(6)!r; end if;
        if posdef and r lt bestR then bestR := r; best := ts; chosen := embset[j]; end if;
        if posdef and r lt basin_thresh then break; end if;  // basin found -> stop searching
    end for;
    error if #best eq 0, "x_to_tau_fast: no valid tau found; increase trials or check x/K";
    if bestR ge basin_thresh then
        printf "WARNING: best search residual %o exceeds basin_thresh %o; polish may not converge\n",
               RealField(8)!bestR, RealField(8)!basin_thresh;
    end if;

    // ---- Phase B: polish the winner at full precision (using the winning embedding) ----
    embH := ComplexEmbedding(K, chosen, CC);
    qsH, dqsH := CoerceQuintics(qs, K, embH, CC);
    // Reconstruct the start in CC from real/imag parts (robust across prec levels).
    tau0_hi := [ CC!Real(t) + CC.1*(CC!Imaginary(t)) : t in best ];
    zs_hi := [ [CC| (Random(-50,50)+CC.1*Random(10,40))/100,
                    (Random(-50,50)+CC.1*Random(10,40))/100 ] : kk in [1..6] ];
    // coarse:=false -> full N from the start (we are in-basin and want full accuracy);
    // lm_floor:=0 -> never abandon the converged solution;
    // tol_pow:=prec-10 -> converge to nearly full precision so recognition succeeds.
    tau_hi, rfinal := FindTauFast(qsH, dqsH, tau0_hi, zs_hi, CC : iters := polish_iters,
                                  verbose := verbose, coarse := false, lm_floor := 0,
                                  tol_pow := prec - 10);
    verify_tau(qsH, tau_hi, CC);             // independent sanity check at full precision
    return SymmetricMatrix(tau_hi), rfinal, embH;  // tau, residual, chosen embedding
end function;
