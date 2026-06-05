// Reference implementation of the inverse map x -> tau (period matrix).
// Strategy: a point z on A_tau maps under the theta embedding to a point that
// must satisfy every quintic q defining the surface A_x. We treat q(Theta(z))
// as a residual depending on tau = [[a,b],[b,c]] and solve by Gauss-Newton.
// InversionFast.m is the optimized, production version of this.
import "HMSurface.m": HMSurface;

prec := 100;
CC := ComplexField(prec);
PI := Pi(CC); ii := CC.1;
P5<x0,x1,x2,x3,x4> := PolynomialRing(CC, 5);
// qs := [q1,q2,q3];  // your three quintics, coefficients coerced K -> CC
//                     // via a FIXED real embedding sqrt(2) |-> 1.41421356...

// (1,5) theta embedding: returns the P^4 point of z in A_tau
function ThetaPt(z, a, b, c, N, CC)
    PI := Pi(CC); ii := CC.1;
    pt := [CC| 0,0,0,0,0];
    for j in [0..4] do
        s := CC!0;
        for n1 in [-N..N], n2 in [-N..N] do
            v1 := CC!n1; v2 := n2 + CC!j/5;
            Q := a*v1^2 + 2*b*v1*v2 + c*v2^2;     // v^t tau v
            L := v1*z[1] + v2*z[2];               // v^t z
            s +:= Exp(PI*ii*Q + 2*PI*ii*L);
        end for;
        pt[j+1] := s;
    end for;
    return pt;
end function;

// scale-clean evaluation: normalize the projective point to unit norm before
// plugging into q, so the residual scale is independent of the theta magnitude.
function QEval(q, pt)
    nrm := Sqrt(&+[ Abs(pt[k])^2 : k in [1..5] ]);
    return Evaluate(q, [ pt[k]/nrm : k in [1..5] ]);
end function;

// Stacked residual vector: q(Theta(z)) for every sample point z and quintic q.
// Vanishes exactly when tau = [[a,b],[b,c]] is the correct period matrix.
function Residual(qs, zs, a, b, c, N, CC)
    R := [CC| ];
    for z in zs do
        pt := ThetaPt(z, a, b, c, N, CC);
        for q in qs do Append(~R, QEval(q, pt)); end for;
    end for;
    return R;
end function;

// Gauss-Newton on the (overdetermined, holomorphic) residual.
// Jacobian is built by central finite differences in a, b, c (step h);
// step solves the normal equations (J* J) del = -J* R.
function FindTau(qs, tau0, zs, N, CC : iters := 40)
    prec := Precision(CC);
    a := tau0[1]; b := tau0[2]; c := tau0[3];
    h := CC!10^(-(prec div 3));               // finite-difference step
    for it in [1..iters] do
        R   := Residual(qs,zs,a,b,c,N,CC);  m := #R;
        Ra  := Residual(qs,zs,a+h,b,c,N,CC); Ra2 := Residual(qs,zs,a-h,b,c,N,CC);
        Rb  := Residual(qs,zs,a,b+h,c,N,CC); Rb2 := Residual(qs,zs,a,b-h,c,N,CC);
        Rc  := Residual(qs,zs,a,b,c+h,N,CC); Rc2 := Residual(qs,zs,a,b,c-h,N,CC);
        J := ZeroMatrix(CC, m, 3);
        for r in [1..m] do
            J[r,1]:=(Ra[r]-Ra2[r])/(2*h);
            J[r,2]:=(Rb[r]-Rb2[r])/(2*h);
            J[r,3]:=(Rc[r]-Rc2[r])/(2*h);
        end for;
        Jstar := ZeroMatrix(CC,3,m);                 // conjugate transpose J*
        for r in [1..m], k in [1..3] do Jstar[k,r] := Conjugate(J[r,k]); end for;
        Rcol := Matrix(CC, m, 1, R);
        // printf "Rank(Jstar*J) = %o\n", Rank(Jstar*J);
        // If the normal matrix is rank-deficient we are at a singular point; stop.
        if Rank(Jstar*J) lt 3 then break; end if;
        del  := -(Jstar*J)^(-1) * (Jstar*Rcol);     // normal equations
        a +:= del[1,1]; b +:= del[2,1]; c +:= del[3,1];
        nrm := Sqrt(&+[Abs(x)^2 : x in R]);
        printf "it %o: |R|=%o  Im tau diag=(%o,%o)\n",
               it, RealField(8)!nrm, RealField(8)!Imaginary(a),
               RealField(8)!Imaginary(c);
        if nrm lt 10^(-(prec div 2)) then break; end if;
    end for;
    return [a,b,c];
end function;

// Independent verification: check that the recovered tau really annihilates the
// quintics on FRESH random points (not used in the fit) and on the theta-null
// (image of the origin, a canonical point of the surface). Asserts on failure.
procedure verify_tau(qs, tau, CC)
    prec := Precision(CC);
    ii := CC.1;
    // fresh points, independent of the Newton run.
    // Nchk must scale with precision (and 1/lm) like the solver's N_need -- a fixed
    // small cap leaves the theta truncation error above the assert tolerance, so a
    // correct tau would fail the check on truncation noise alone.
    lm := ((Imaginary(tau[1])+Imaginary(tau[3]))
           - Sqrt((Imaginary(tau[1])-Imaginary(tau[3]))^2+4*Imaginary(tau[2])^2))/2;
    Nchk := Ceiling(Sqrt(prec*Log(10.0)/(Real(Pi(CC))*Max(lm, Real(CC!0.005)))))+3;
    zs2 := [ [CC| (Random(-70,70)+ii*Random(15,45))/100,
                (Random(-70,70)+ii*Random(15,45))/100 ] : k in [1..10] ];
    res := Residual(qs, zs2, tau[1], tau[2], tau[3], Nchk, CC);
    printf "max |q(P)| on fresh points: %o\n", Max([Abs(x): x in res]);
    // the theta-null point (image of the origin) is a canonical point of the surface:
    res0 := [ QEval(q, ThetaPt([CC|0,0], tau[1],tau[2],tau[3], Nchk, CC)) : q in qs ];
    printf "max |q(theta-null)|: %o\n", Max([Abs(x): x in res0]);
    assert Max([Abs(x): x in res]) lt 10^(-(prec div 3));
    assert Max([Abs(x): x in res0]) lt 10^(-(prec div 3));
    return;
end procedure;

// !! TODO = make all the constants into parameters
// x in P^3 to tau such that A_tau = A_x
function x_to_tau(x)
    k := Rationals();                      // GF(p) with p prime is much faster;
    P := ProjectiveSpace(k, 4);
    F := HorrocksMumfordBundle(P);
    A := HMSurface(F, x);
    basis := MinimalBasis(Ideal(A));
    qs := [b : b in basis | Degree(b) eq 5];      // the 3 defining quintics
    assert #qs eq 3;
    // CC := ComplexField(200);
    // Sample points on A_tau, drawn from a reduced box in the upper half plane.
    zs := [ [CC| (Random(-50,50)+ii*Random(10,40))/100,
             (Random(-50,50)+ii*Random(10,40))/100 ] : k in [1..6] ];
    best := []; bestR := Infinity();
    // Random restarts: many initial tau0 guesses, keep the best pos-def solution.
    for trial in [1..20] do
        tau0 := [ (Random(-50,50)+ii*Random(20,80))/100,    // a, in a reduced box
                (Random(-25,25))/100,                       // b small
                (Random(-50,50)+ii*Random(20,80))/100 ];    // c
        t := FindTau(qs, tau0, zs, 12, CC);
        // Keep only valid period matrices: Im(tau) positive definite.
        posdef := Imaginary(t[1]) gt 0 and
                Imaginary(t[1])*Imaginary(t[3]) - Imaginary(t[2])^2 gt 0;
        r := Sqrt(&+[Abs(x)^2 : x in Residual(qs,zs,t[1],t[2],t[3],12,CC)]);
        if posdef and r lt bestR then bestR := r; best := t; end if;
    end for;
    verify_tau(qs, best, CC);
    return SymmetricMatrix(best), bestR;     // [a,b,c] -> [[a,b],[b,c]]
end function;
