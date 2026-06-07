// Fast inverse for the (1,11) pipeline: pencil <v,w> in V_+  ->  period tau, with
// A_tau = the (1,11) abelian surface cut out by the 22 quadrics of the pencil. This is
// the (1,11) analogue of InversionFast.m (analytic Jacobian + Levenberg-Marquardt +
// adaptive truncation + two-phase precision). The convention (ThetaPt11 vs the R5
// quadrics) is verified consistent in GP11.m, so a random tau in a basin zeroes the
// residual q(Theta_tau(z)) over the surface's quadrics.
import "GP11.m": IntertwiningMatrix, PencilFromTorsion;

D := 11;

// Theta values + the three tau-derivatives (d/da, d/db, d/dc) for all 11 (1,11) theta
// coordinates, in ONE lattice pass (each term's tau-derivative is pi*i times a quadratic
// monomial in v times the term). Direct analogue of ThetaJet for (1,5).
function ThetaJet11(z, a, b, c, N, CC)
    PI := Pi(CC); ii := CC.1;
    th := [CC| 0 : i in [1..D]]; tha := th; thb := th; thc := th;
    for j in [0..D-1] do
        sj := CC!0; sa := CC!0; sb := CC!0; sc := CC!0;
        for n1 in [-N..N], n2 in [-N..N] do
            v1 := CC!n1; v2 := n2 + CC!j/D;
            E := Exp(PI*ii*(a*v1^2 + 2*b*v1*v2 + c*v2^2) + 2*PI*ii*(v1*z[1] + v2*z[2]));
            sj +:= E; sa +:= (PI*ii*v1^2)*E; sb +:= (PI*ii*2*v1*v2)*E; sc +:= (PI*ii*v2^2)*E;
        end for;
        th[j+1] := sj; tha[j+1] := sa; thb[j+1] := sb; thc[j+1] := sc;
    end for;
    return th, tha, thb, thc;
end function;

// The 22 quadrics of the pencil <v,w> (length-6 over CC) as polynomials in P11 = CC[X0..X10],
// plus their 11 partials (precomputed once). Q^u_j = sum_{i=0}^5 u_i X_{j+i} X_{j-i}.
function QuadricsCC(v, w, P11)
    X := [P11.i : i in [1..D]];
    qs := [ &+[ u[i+1] * X[((j+i) mod D)+1] * X[((j-i) mod D)+1] : i in [0..5] ]
            : j in [0..D-1], u in [v, w] ];
    return qs, [ [ Derivative(q, l) : l in [1..D] ] : q in qs ];
end function;

// Residual vector R and analytic Jacobian J (w.r.t. a,b,c). For each z, dehomogenize the
// theta point at its largest coordinate (a safe chart, fixed within this evaluation, so the
// quotient-rule derivatives are exact), then chain-rule the quadric gradients. Mirrors ResJac.
function ResJac11(qs, dqs, zs, a, b, c, N, CC)
    Rseq := [CC| ]; JA := [CC| ]; JB := [CC| ]; JC := [CC| ];
    for z in zs do
        th, tha, thb, thc := ThetaJet11(z, a, b, c, N, CC);
        _, k := Max([Abs(th[l]) : l in [1..D]]);
        mu := th[k]; ma := tha[k]; mb := thb[k]; mc := thc[k];
        Pt := [ th[l]/mu : l in [1..D] ];
        Pa := [ (tha[l]*mu - th[l]*ma)/mu^2 : l in [1..D] ];
        Pb := [ (thb[l]*mu - th[l]*mb)/mu^2 : l in [1..D] ];
        Pc := [ (thc[l]*mu - th[l]*mc)/mu^2 : l in [1..D] ];
        for i in [1..#qs] do
            Append(~Rseq, Evaluate(qs[i], Pt));
            grad := [ Evaluate(dqs[i][l], Pt) : l in [1..D] ];
            Append(~JA, &+[ grad[l]*Pa[l] : l in [1..D] ]);
            Append(~JB, &+[ grad[l]*Pb[l] : l in [1..D] ]);
            Append(~JC, &+[ grad[l]*Pc[l] : l in [1..D] ]);
        end for;
    end for;
    m := #Rseq; R := Matrix(CC, m, 1, Rseq); J := ZeroMatrix(CC, m, 3);
    for r in [1..m] do J[r,1] := JA[r]; J[r,2] := JB[r]; J[r,3] := JC[r]; end for;
    return R, J;
end function;

function ResNorm11(qs, dqs, zs, a, b, c, N, CC)
    R := ResJac11(qs, dqs, zs, a, b, c, N, CC);
    return Sqrt(&+[ Abs(R[i,1])^2 : i in [1..Nrows(R)] ]);
end function;

// Levenberg-Marquardt solve for tau from tau0 (identical control logic to FindTauFast).
function FindTauFast11(qs, dqs, tau0, zs, CC : iters := 60, verbose := false,
                       stag_floor := 1e-3, stag_rel := 1e-3, stag_win := 4, coarse := true,
                       lm_floor := 1e-2, tol_pow := 0)
    prec := Precision(CC);
    a := tau0[1]; b := tau0[2]; c := tau0[3];
    lambda := CC!1e-4; nrm := 10^9;
    tol := tol_pow gt 0 select Real(CC!10^(-tol_pow)) else Real(CC!10^(-Floor(3*prec/7)));
    prevnrm := 10^18; stag := 0;
    for it in [1..iters] do
        lm := ((Imaginary(a)+Imaginary(c)) - Sqrt((Imaginary(a)-Imaginary(c))^2+4*Imaginary(b)^2))/2;
        if lm_floor gt 0 and lm lt lm_floor then
            if verbose then printf "  [Im(tau) eigenvalue %o below floor, abandoning]\n", RealField(6)!lm; end if; break;
        end if;
        N_need := Ceiling(Sqrt(prec*Log(10.0)/(Real(Pi(CC))*Max(lm, Real(CC!0.005)))))+3;
        if not coarse then N := N_need;
        elif nrm gt 10^(-2) then N := Min(N_need, Max(8, N_need div 3));
        elif nrm gt 10^(-6) then N := Min(N_need, Max(12, (2*N_need) div 3));
        else N := N_need; end if;
        R, J := ResJac11(qs, dqs, zs, a, b, c, N, CC); m := Nrows(R);
        nrm := Sqrt(&+[ Abs(R[i,1])^2 : i in [1..m] ]);
        if verbose then printf "it %o: N=%o lam_min=%o |R|=%o\n", it, N, RealField(6)!lm, RealField(8)!nrm; end if;
        if nrm lt tol then break; end if;
        if nrm gt stag_floor and nrm ge prevnrm*(1 - stag_rel) then stag +:= 1; else stag := 0; end if;
        if stag ge stag_win then
            if verbose then printf "  [stalled at |R|=%o, abandoning]\n", RealField(8)!nrm; end if; break;
        end if;
        prevnrm := nrm;
        Js := ZeroMatrix(CC, 3, m);
        for i in [1..m], k in [1..3] do Js[k,i] := Conjugate(J[i,k]); end for;
        JsJ := Js*J;
        if Rank(JsJ) lt 3 then
            if verbose then printf "  [singular normal matrix, abandoning]\n"; end if; break;
        end if;
        accepted := false;
        step := -(JsJ)^(-1)*(Js*R);
        an := a + step[1,1]; bn := b + step[2,1]; cn := c + step[3,1];
        if Imaginary(an) gt 0 and Imaginary(an)*Imaginary(cn) - Imaginary(bn)^2 gt 0 then
            nn := ResNorm11(qs, dqs, zs, an, bn, cn, N, CC);
            if nn lt nrm then a := an; b := bn; c := cn; lambda := CC!1e-4; accepted := true; end if;
        end if;
        Dm := DiagonalMatrix([ JsJ[i,i] : i in [1..3] ]);
        while not accepted and Abs(lambda) lt 1e18 do
            step := -(JsJ + lambda*Dm)^(-1)*(Js*R);
            an := a + step[1,1]; bn := b + step[2,1]; cn := c + step[3,1];
            if Imaginary(an) gt 0 and Imaginary(an)*Imaginary(cn) - Imaginary(bn)^2 gt 0 then
                nn := ResNorm11(qs, dqs, zs, an, bn, cn, N, CC);
                if nn lt nrm then a := an; b := bn; c := cn; lambda := lambda/4; accepted := true;
                else lambda := lambda*3; end if;
            else lambda := lambda*3; end if;
        end while;
    end for;
    return [a,b,c], nrm;
end function;

// Two-phase inversion: locate the basin with cheap reduced-precision restarts, then polish
// the winner at full precision. Input: the pencil v,w (length-6 over CC). Returns <tau, |R|>.
function InvertGP11Fast(v, w, CC : trials := 100, search_prec := 0, verbose := false,
                        basin_thresh := 1e-6, polish_iters := 40)
    prec := Precision(CC);
    sp := search_prec gt 0 select search_prec else Max(30, Min(prec, 40));
    CCs := ComplexFieldExtra(sp); iis := CCs.1;
    P11s := PolynomialRing(CCs, D);
    qsS, dqsS := QuadricsCC([CCs!c : c in v], [CCs!c : c in w], P11s);
    zs := [ [CCs| (Random(-50,50)+iis*Random(10,40))/100, (Random(-50,50)+iis*Random(10,40))/100 ] : kk in [1..6] ];
    Nscore := Max(14, Ceiling(Sqrt(sp/4)));
    best := []; bestR := Infinity();
    for trial in [1..trials] do
        tau0 := [ (Random(-200,200)+iis*Random(10,120))/100,
                  (Random(-150,150)+iis*Random(-50,50))/100,
                  (Random(-250,250)+iis*Random(10,250))/100 ];
        ts, _ := FindTauFast11(qsS, dqsS, tau0, zs, CCs : iters := 50, verbose := verbose);
        pd := Imaginary(ts[1]) gt 0 and Imaginary(ts[1])*Imaginary(ts[3]) - Imaginary(ts[2])^2 gt 0;
        r := pd select ResNorm11(qsS, dqsS, zs, ts[1], ts[2], ts[3], Nscore, CCs) else Infinity();
        if pd and r lt bestR then bestR := r; best := ts; end if;
        if verbose then printf "trial %o: posdef=%o |R|=%o\n", trial, pd, pd select RealField(6)!r else r; end if;
        if pd and r lt basin_thresh then break; end if;
    end for;
    error if #best eq 0, "InvertGP11Fast: no positive-definite tau found; increase trials";
    if bestR ge basin_thresh then
        printf "WARNING: best search residual %o exceeds basin_thresh %o\n", RealField(8)!bestR, RealField(8)!basin_thresh;
    end if;
    P11 := PolynomialRing(CC, D);
    qsH, dqsH := QuadricsCC(v, w, P11);
    tau0_hi := [ CC!Real(t) + CC.1*(CC!Imaginary(t)) : t in best ];
    zs_hi := [ [CC| (Random(-50,50)+CC.1*Random(10,40))/100, (Random(-50,50)+CC.1*Random(10,40))/100 ] : kk in [1..6] ];
    tau_hi, rfinal := FindTauFast11(qsH, dqsH, tau0_hi, zs_hi, CC : iters := polish_iters,
                                    verbose := verbose, coarse := false, lm_floor := 0, tol_pow := prec - 10);
    return SymmetricMatrix(tau_hi), rfinal;
end function;
