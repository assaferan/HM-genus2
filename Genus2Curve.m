// Top-level driver: recover the genus-2 curve C with Jac(C) ~ A_x (a rational
// cyclic 5-isogeny), for x in P^3. CHIMP (sibling dir ../CHIMP) supplies
// ComplexFieldExtra / ReduceSmallPeriodMatrix / ThetaDerivatives.
AttachSpec("../CHIMP/CHIMP.spec");
import "InversionFast.m": x_to_tau_fast;

// Recognize a complex number as a rational (degree-1 algebraic number).
// Returns the rational and the (tiny) residual |a*z - b| for sanity checking.
function RecognizeRational(z)
    mp := MinimalPolynomial(z, 1);
    q := -Coefficient(mp, 0) / Coefficient(mp, 1);
    res := Abs(Evaluate(ChangeRing(mp, Parent(z)), z));
    return q, res;
end function;

// Period matrix tau (of the principally polarized [tau|I]) -> genus-2 curve over Q,
// reconstructed from its rational ABSOLUTE Igusa invariants. This replicates the
// numeric core of CHIMP's AlgebraizedInvariantsG2 (theta derivatives at the six odd
// two-torsion points -> branch points -> sextic -> Igusa invariants), then recognizes
// the invariants over Q. We avoid CHIMP's ReconstructCurve here because it algebraizes
// the *twisted sextic coefficients*, whose height dwarfs that of the absolute
// invariants and overflows LLL's height bound at usable precision. The curve is thus
// determined up to quadratic twist -- harmless here, since a rational 5-isogeny is
// preserved by quadratic twist (the stable line in J[5] survives tensoring by a
// quadratic character).
function Genus2CurveFromTau(tau)
    CC := BaseRing(tau);
    P := HorizontalJoin(tau, IdentityMatrix(CC, 2));
    taunew, gamma := ReduceSmallPeriodMatrix(tau);     // Siegel-reduce for stable theta
    Am := Submatrix(gamma,1,1,2,2); Bm := Submatrix(gamma,1,3,2,2);
    Cm := Submatrix(gamma,3,1,2,2); Dm := Submatrix(gamma,3,3,2,2);
    Pnew := P * Transpose(BlockMatrix([[Am,Bm],[Cm,Dm]]));
    P2inew := Submatrix(Pnew,1,3,2,2)^(-1);
    // the six odd two-torsion points (1/2)(taunew*a + b)
    half := func< a, b | (1/2)*taunew*Transpose(Matrix(CC,[a])) + (1/2)*Transpose(Matrix(CC,[b])) >;
    ws := [ half([0,1],[0,1]), half([0,1],[1,1]), half([1,0],[1,0]),
            half([1,0],[1,1]), half([1,1],[0,1]), half([1,1],[1,0]) ];
    theta_derss := [ ThetaDerivatives(taunew, w) : w in ws ];
    Hs := [ Matrix(CC,[td]) * P2inew : td in theta_derss ];
    rats := [];
    for H in Hs do
        seq := Eltseq(H); add := true;
        if Abs(seq[2]) lt Abs(seq[1]) and Abs(seq[2]/seq[1])^2 lt CC`epscomp then add := false; end if;
        if add then Append(~rats, -seq[1]/seq[2]); end if;
    end for;
    error if #rats ne 6, "Genus2CurveFromTau: expected 6 branch points, got", #rats;
    RCC := PolynomialRing(CC);
    fCC := &*[ RCC.1 - r : r in rats ];
    ICCn := WPSNormalizeCC([2,4,6,8,10], IgusaInvariants(fCC));   // absolute invariants
    QI := [Rationals()| ];
    for z in ICCn do
        q, res := RecognizeRational(z);
        error if res gt 1e-20,
            "Genus2CurveFromTau: Igusa invariant not rational (curve over a number field?); increase Prec";
        Append(~QI, q);
    end for;
    return HyperellipticCurveFromIgusaInvariants(QI);
end function;

function Genus2Curve(x : Prec := 300)
    CC := ComplexFieldExtra(Prec);
    tau := x_to_tau_fast(x, CC);            // numeric inverse: x -> period matrix tau
    return Genus2CurveFromTau(tau);         // tau -> curve over Q via rational invariants
end function;
