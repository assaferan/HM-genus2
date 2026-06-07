// (1,11) analogue of the Horrocks-Mumford pipeline, from Gross-Popescu, "The moduli
// space of (1,11)-polarized abelian surfaces is unirational" (arXiv math/9902017; the
// (1,11) details that [GP1] "Equations of (1,d)-polarized abelian surfaces", Math. Ann.
// 310 (1998), works out in general).  Moduli A_11^lev is birational to the KLEIN CUBIC
//   K = { x0^2 x1 + x1^2 x2 + x2^2 x3 + x3^2 x4 + x4^2 x0 = 0 } subset P^4
// (the unique PSL2(Z/11)-invariant cubic).  Target pipeline, paralleling (1,5):
//   moduli point  -->  (1,11) abelian surface A subset P^10  -->  period tau
//   -->  principally polarized [tau|I] (11-isogenous to A)  -->  genus-2 curve C with a
//   rational 11-isogeny.
//
// EXPLICIT CONSTRUCTION (GP2):
//  * Heisenberg H_11 acts on P^10 = P(C[Z/11]) by  sigma(x_i)=x_{i-1},
//    tau11(x_i)=xi^{-i} x_i  (xi = exp(2 pi i/11)),  iota(x_i)=x_{-i}.
//  * S^2(V) (the 66 quadrics) = 6 copies of an 11-dim H_11-rep, organized by the 6x11
//    matrix  (R5)_{ij} = x_{j+i} x_{j-i}  (0<=i<=5, 0<=j<=10, indices mod 11).  A vector
//    v=(v_0..v_5) in V_+ = C^6 gives the 11 quadrics  Q^v_j = sum_i v_i x_{j+i} x_{j-i}.
//    A (1,11) surface A is cut out by the 22 quadrics of a PENCIL <v,w> in V_+, i.e. a
//    point of Gr(2,6) (image of A under Theta_11; im(Theta) is birational to K).
//  * Klein cubic = Pfaffian of the 6x6 skew matrix M (linear in x_0..x_4) below.
//  * The intertwining operator is the 6x6 skew matrix S (quadratic in x_1..x_5); the
//    sextic D_2 = {Pf(S)=0} subset P^4 is the locus of odd 2-torsion points, and the
//    EXPLICIT map  P |-> ker S(P)  sends P in D_2 to the pencil <v,w> of A.  (The Klein
//    cubic <-> Gr(2,6) identification is the *indirect* Fano-Iskovskih birational map.)
//
// STATUS: the algebraic construction (Klein cubic, M, S, R5 quadrics, surface from a
// 2-torsion point) is implemented + self-tested below.  TODO: port FindTauFast to invert
// the (1,11) theta map against these quadrics, then reuse Genus2Curve.m / Reduction.m
// verbatim (with "5" -> "11": good reduction at 11, L_p mod 11, rational 11-isogeny).

D := 11;

// ---- numeric (1,11) theta embedding: A_tau -> P^10, analogue of (1,5) ThetaPt --------
function ThetaPt11(z, tau, N, CC)
    PI := Pi(CC); ii := CC.1; a := tau[1][1]; b := tau[1][2]; c := tau[2][2];
    pt := [CC| 0 : j in [1..D]];
    for j in [0..D-1] do
        s := CC!0;
        for n1 in [-N..N], n2 in [-N..N] do
            v1 := CC!n1; v2 := n2 + CC!j/D;
            s +:= Exp(PI*ii*(a*v1^2 + 2*b*v1*v2 + c*v2^2) + 2*PI*ii*(v1*z[1] + v2*z[2]));
        end for;
        pt[j+1] := s;
    end for;
    return pt;
end function;

// ---- Heisenberg H_11 on P^10 (the symmetry the (1,11) equations respect) --------------
function HeisenbergGenerators11(CC)
    xi := Exp(2*Pi(CC)*CC.1/D);
    sigma := Matrix(CC, D, D, [ <i, ((i-2) mod D)+1, CC!1> : i in [1..D] ]);  // x_i -> x_{i-1}
    tau11 := DiagonalMatrix(CC, [ xi^(-(i-1)) : i in [1..D] ]);               // x_i -> xi^{-i} x_i
    iota  := Matrix(CC, D, D, [ <i, ((-(i-1)) mod D)+1, CC!1> : i in [1..D] ]);// x_i -> x_{-i}
    return sigma, tau11, iota;
end function;

// ---- Klein cubic threefold and its Pfaffian matrix M (P^4, coords x0..x4) -------------
// M is 6x6 skew, linear in x0..x4, with Pf(M) = x0^2 x1 + ... + x4^2 x0 (the Klein cubic).
function KleinCubicMatrix(x)   // x = [x0,x1,x2,x3,x4]
    x0,x1,x2,x3,x4 := Explode(x); R := Universe(x); z := R!0;
    return Matrix(R, 6, 6, [
        [  z,  x0,  x2,  x1,  x4,  x3 ],
        [-x0,   z,  x4,   z,   z, -x2 ],
        [-x2, -x4,   z,   z,  x1,   z ],
        [-x1,   z,   z,   z, -x3,  x0 ],
        [-x4,   z, -x1,  x3,   z,   z ],
        [-x3,  x2,   z, -x0,   z,   z ] ]);
end function;

// ---- Intertwining operator S (P^4, coords x1..x5) ; D_2 = {Pf(S)=0} (odd 2-torsion) ---
function IntertwiningMatrix(x)  // x = [x1,x2,x3,x4,x5]
    x1,x2,x3,x4,x5 := Explode(x); R := Universe(x); z := R!0;
    return Matrix(R, 6, 6, [
        [    z,  x1^2,  x2^2,  x3^2,  x4^2,   x5^2 ],
        [-x1^2,     z, x1*x3, x2*x4, x3*x5, -x4*x5 ],
        [-x2^2,-x1*x3,     z, x1*x5,-x2*x5, -x3*x4 ],
        [-x3^2,-x2*x4,-x1*x5,     z,-x1*x4, -x2*x3 ],
        [-x4^2,-x3*x5, x2*x5, x1*x4,     z, -x1*x2 ],
        [-x5^2, x4*x5, x3*x4, x2*x3, x1*x2,      z ] ]);
end function;

// ---- the 11 quadrics in P^10 (coords X0..X10) attached to v=(v0..v5) in V_+ -----------
// Q^v_j = sum_{i=0}^5 v_i X_{j+i} X_{j-i}, indices mod 11 (the rows of (R5)_{ij}=X_{j+i}X_{j-i}).
function QuadricsFromVector(v, X)   // v: length 6 (over the ring of X); X: [X0..X10]
    R := Universe(X);
    Qs := [R| ];
    for j in [0..D-1] do
        Append(~Qs, &+[ v[i+1] * X[((j+i) mod D)+1] * X[((j-i) mod D)+1] : i in [0..5] ]);
    end for;
    return Qs;
end function;

// 22 quadrics of the surface from a pencil <v,w> in V_+ (the two basis vectors).
function QuadricsFromPencil(v, w, X)
    return QuadricsFromVector(v, X) cat QuadricsFromVector(w, X);
end function;

// The 22 quadrics in P^10 (polynomial ring R10, coords X0..X10) cutting out the (1,11)
// surface whose odd 2-torsion point is Ptors in D_2: the pencil <v,w> = ker S(Ptors) in
// V_+ (Theta: P |-> ker S(P)), then Q^v_j, Q^w_j via the R5 matrix.  Ptors is a 5-vector
// (x1..x5) over the base field of R10; it must lie on D_2 \ D_1 (so ker S has dim 2).
function SurfaceQuadrics(Ptors, R10)
    F := BaseRing(R10);
    S := IntertwiningMatrix([F!c : c in Ptors]);
    K := Kernel(S);
    error if Dimension(K) ne 2, "Ptors not in D_2 \\ D_1 (ker S has dim", Dimension(K), "!= 2)";
    bw := Basis(K);
    v := [R10!c : c in Eltseq(bw[1])];
    w := [R10!c : c in Eltseq(bw[2])];
    return QuadricsFromPencil(v, w, [R10.i : i in [1..11]]);
end function;

// Pencil <v,w> = ker S(Ptors) over an exact field (for a 2-torsion point Ptors on D_2).
// Returned as two 6-vectors; embed into CC for the theta inversion.
function PencilFromTorsion(Ptors)
    S := IntertwiningMatrix(Ptors);
    K := Kernel(S);
    error if Dimension(K) ne 2, "Ptors not on D_2 \\ D_1 (ker S has dim", Dimension(K), ")";
    bw := Basis(K);
    return Eltseq(bw[1]), Eltseq(bw[2]);
end function;

// ---- (1,11) theta inversion: find tau with Theta_tau(A) on the surface's quadrics -----
// The surface A is given by a pencil v,w in V_+ (length-6 vectors over CC). A point z lies
// on A_tau, so its (1,11) theta image must satisfy all 22 quadrics Q^v_j, Q^w_j. We solve
// for tau = [[a,b],[b,c]] by (complex) Gauss-Newton on these residuals over sample z's.
// Mirrors Inversion.m for (1,5); the residual is holomorphic in tau (theta point is
// normalized by a FIXED coordinate, not by abs), so the Jacobian can be taken over C.

function ThetaMatrix11(abc)              // abc = [a,b,c] -> symmetric 2x2
    return Matrix(Parent(abc[1]), 2, 2, [abc[1], abc[2], abc[2], abc[3]]);
end function;

function ResidualAtZ11(v, w, z, abc, N, CC)
    P := ThetaPt11(z, ThetaMatrix11(abc), N, CC);  // raw theta point (residual stays holomorphic in tau)
    res := [CC| ];
    for u in [v, w] do
        for j in [0..D-1] do
            Append(~res, &+[ u[i+1] * P[((j+i) mod D)+1] * P[((j-i) mod D)+1] : i in [0..5] ]);
        end for;
    end for;
    return res;
end function;

function Residual11(v, w, zs, abc, N, CC)
    return &cat[ ResidualAtZ11(v, w, z, abc, N, CC) : z in zs ];
end function;

// Complex Gauss-Newton with Levenberg-Marquardt damping. Returns <abc, residual norm>.
function FindTau11(v, w, zs, abc0, N, CC : iters := 80, lambda0 := 1e-3, verbose := false)
    RR := RealField(Precision(CC));
    abc := abc0;
    r := Residual11(v, w, zs, abc, N, CC);
    nr := Sqrt(&+[ Abs(t)^2 : t in r ]);
    lambda := RR!lambda0;
    eps := (RR!10)^(-(Precision(CC) div 3));
    for it in [1..iters] do
        cols := [];
        for k in [1..3] do
            abck := abc; abck[k] +:= eps;
            rk := Residual11(v, w, zs, abck, N, CC);
            Append(~cols, [ (rk[l] - r[l]) / eps : l in [1..#r] ]);
        end for;
        J := Matrix(CC, #r, 3, [ cols[k][l] : l in [1..#r], k in [1..3] ]);
        Jh := Transpose(J); Jh := Matrix(CC, 3, #r, [ ComplexConjugate(Jh[i][j]) : i in [1..3], j in [1..#r] ]);
        H := Jh * J; g := Jh * Matrix(CC, #r, 1, r);
        ok := false; step := [CC|0,0,0];
        for tries in [1..12] do
            Hd := H + lambda * IdentityMatrix(CC, 3);
            d := -(Hd^(-1)) * g;
            cand := [ abc[k] + d[k][1] : k in [1..3] ];
            rc := Residual11(v, w, zs, cand, N, CC);
            nc := Sqrt(&+[ Abs(t)^2 : t in rc ]);
            if nc lt nr then abc := cand; r := rc; nr := nc; lambda /:= 2; ok := true; break;
            else lambda *:= 3; end if;
        end for;
        if verbose then printf "   it %o: |r|=%o lambda=%o\n", it, ChangePrecision(nr,6), ChangePrecision(lambda,3); end if;
        if not ok or nr lt eps*100 then break; end if;
    end for;
    return abc, nr;
end function;

// Invert with random restarts; keep the best tau with positive-definite imaginary part.
function InvertGP11(v, w, CC : trials := 40, N := 8, nz := 6, verbose := false)
    RR := RealField(Precision(CC)); SetSeed(1);
    zs := [ [ (CC!Random(-500,500) + CC.1*CC!Random(-500,500))/1000 : i in [1..2] ] : k in [1..nz] ];
    best := [CC|0,0,0]; bestnr := (RR!10)^9; got := false;
    for t in [1..trials] do
        a := CC!Random(-300,300)/1000 + CC.1*(CC!Random(200,1200)/1000);
        b := CC!Random(-300,300)/1000 + CC.1*CC!Random(-200,200)/1000;
        c := CC!Random(-300,300)/1000 + CC.1*(CC!Random(200,1200)/1000);
        abc, nr := FindTau11(v, w, zs, [a,b,c], N, CC);
        i11 := Imaginary(abc[1]); i12 := Imaginary(abc[2]); i22 := Imaginary(abc[3]);
        posdef := i11 gt 0 and i11*i22 - i12^2 gt 0;
        if verbose then printf "  trial %o: |r|=%o posdef=%o\n", t, ChangePrecision(nr,4), posdef; end if;
        if nr lt bestnr and posdef then best := abc; bestnr := nr; got := true; end if;
    end for;
    return best, bestnr, got;
end function;

// im(Theta_11) subset Gr(2,6): the 5 linear Plucker relations cutting the moduli 3-fold
// out of Gr(2,6) (Theorem 2.6).  p_ij are Plucker coords of the pencil <v,w> in V_+.
//   p23=-p15, p26=p13, p14=-p35, p16=p45, p46=-p12.
PluckerRelationsX := [
    "p23 + p15", "p26 - p13", "p14 + p35", "p16 - p45", "p46 + p12" ];

// ---- self-test of the algebra --------------------------------------------------------
procedure GP11SelfTest()
    R4<x0,x1,x2,x3,x4> := PolynomialRing(Rationals(), 5);
    M := KleinCubicMatrix([x0,x1,x2,x3,x4]);
    pf := Pfaffian(M);
    klein := x0^2*x1 + x1^2*x2 + x2^2*x3 + x3^2*x4 + x4^2*x0;
    printf "Pf(M) = Klein cubic : %o\n", pf eq klein or pf eq -klein;
    S := IntertwiningMatrix([x0,x1,x2,x3,x4]);   // (x1..x5) := (x0..x4) here
    printf "S is skew-symmetric : %o\n", Transpose(S) eq -S;
    printf "Pf(S) is a sextic   : %o (degree %o)\n", TotalDegree(Pfaffian(S)) eq 6, TotalDegree(Pfaffian(S));
    // a generic pencil gives 22 quadrics spanning an H_11-invariant space (degree 2 check)
    R10 := PolynomialRing(Rationals(), 11);
    X := [R10.i : i in [1..11]];
    v := [Rationals()| 1,0,0,0,0,0 ]; w := [Rationals()| 0,1,0,0,0,0 ];
    Qs := QuadricsFromPencil(v, w, X);
    printf "pencil -> %o quadrics, all degree 2 : %o\n", #Qs, &and[TotalDegree(q) eq 2 : q in Qs];
end procedure;
