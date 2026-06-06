// Top-level driver: recover the genus-2 curve C with Jac(C) ~ A_x (a rational
// cyclic 5-isogeny) for a point x in P^3(K), where K is Q (the default) or a number
// field. CHIMP (sibling dir ../CHIMP) supplies ComplexFieldExtra / theta machinery.
//
//   Genus2Curve([1,2,3,4]);                                  // over Q
//   K<s2> := QuadraticField(2);
//   Genus2Curve([K| 1,s2,3,4] : K := K);                     // over Q(sqrt2)
//
// Over a number field only some complex embeddings of K yield a positive-definite
// period tau (which ones is x-dependent), so Embedding := 0 (default) auto-tries the
// embeddings and uses the first that works; Embedding := i forces the i-th. Different
// (working) embeddings give Galois-conjugate curves. The curve is reconstructed from
// its absolute Igusa invariants (recognized in K) and, when K != Q, descended to K
// via Mestre's algorithm. A model over K exists iff the Mestre conic has a K-point;
// the Igusa invariants always lie in K (the field of moduli).
AttachSpec("../CHIMP/CHIMP.spec");
import "InversionFast.m": x_to_tau_fast;

// Recognize a complex value as an element of K (degree <= [K:Q]) under embedding emb.
// Returns <element of K, recognition residual>.
function RecognizeInField(z, K, emb)
    degK := Type(K) eq FldRat select 1 else Degree(K);
    mp := MinimalPolynomial(z, degK);
    res := Abs(Evaluate(ChangeRing(mp, Parent(z)), z));
    if Degree(mp) eq 1 then
        return K!(-Coefficient(mp,0)/Coefficient(mp,1)), res;
    end if;
    rts := [ r[1] : r in Roots(mp, K) ];
    error if #rts eq 0, "RecognizeInField: minimal polynomial has no root in K";
    best := rts[1]; bd := Abs(emb(best) - z);
    for r in rts do
        if Abs(emb(r) - z) lt bd then bd := Abs(emb(r) - z); best := r; end if;
    end for;
    return best, res;
end function;

// Numerical absolute (WPS-normalized) Igusa invariants of the ppav [tau|I].
// Replicates the numeric core of CHIMP's AlgebraizedInvariantsG2: theta derivatives
// at the six odd two-torsion points -> branch points -> sextic -> Igusa invariants.
function IgusaNumericFromTau(tau)
    CC := BaseRing(tau);
    P := HorizontalJoin(tau, IdentityMatrix(CC, 2));
    taunew, gamma := ReduceSmallPeriodMatrix(tau);     // Siegel-reduce for stable theta
    Am := Submatrix(gamma,1,1,2,2); Bm := Submatrix(gamma,1,3,2,2);
    Cm := Submatrix(gamma,3,1,2,2); Dm := Submatrix(gamma,3,3,2,2);
    Pnew := P * Transpose(BlockMatrix([[Am,Bm],[Cm,Dm]]));
    P2inew := Submatrix(Pnew,1,3,2,2)^(-1);
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
    error if #rats ne 6, "IgusaNumericFromTau: expected 6 branch points, got", #rats;
    RCC := PolynomialRing(CC);
    fCC := &*[ RCC.1 - r : r in rats ];
    return WPSNormalizeCC([2,4,6,8,10], IgusaInvariants(fCC));
end function;

// Period matrix tau (of the ppav [tau|I]) -> genus-2 curve over K, reconstructed from
// its Igusa invariants recognized in K (embedding emb). When the model lands on an
// extension of K, descend to K via Mestre (Igusa-Clebsch). Reduce when K = Q. The
// curve is determined up to quadratic twist -- harmless for the 5-isogeny, which is
// twist-invariant. We avoid CHIMP's ReconstructCurve because it algebraizes the
// twisted sextic coefficients, whose height dwarfs that of the absolute invariants.
// Absolute Igusa invariants of [tau|I] recognized in K (no curve built). This is the
// cheap part to use when screening many x -- the Mestre descent is only needed to
// produce a model, which we defer to survivors.
function IgusaInvariantsInK(tau, K, emb)
    ICCn := IgusaNumericFromTau(tau);
    QI := [K| ];
    for z in ICCn do
        el, res := RecognizeInField(z, K, emb);
        error if res gt 1e-20,
            "IgusaInvariantsInK: Igusa invariant not recognized in K; increase Prec (or check the embedding)";
        Append(~QI, el);
    end for;
    return QI;
end function;

function Genus2CurveFromTau(tau, K, emb)
    QI := IgusaInvariantsInK(tau, K, emb);
    C := HyperellipticCurveFromIgusaInvariants(QI);
    if not (BaseRing(C) cmpeq K) then          // descend to K via Mestre's algorithm
        IC := [K| c : c in IgusaClebschInvariants(C) ];
        C := HyperellipticCurveFromIgusaClebsch(IC);
    end if;
    if Type(BaseRing(C)) eq FldRat then        // Magma's reduction only works over Q
        try C := ReducedMinimalWeierstrassModel(C); catch e ; end try;
    end if;
    return C, QI;
end function;

// Main driver. K defaults to Q; pass K := <number field> and x over K for the
// number-field case. Embedding picks the complex embedding (number-field case).
function Genus2Curve(x : K := Rationals(), Embedding := 0, Prec := 300,
                     Trials := 100, Verbose := false)
    CC := ComplexFieldExtra(Prec);
    tau, _, emb := x_to_tau_fast(x, CC : K := K, Embedding := Embedding,
                                 trials := Trials, verbose := Verbose);
    return Genus2CurveFromTau(tau, K, emb);   // returns <curve over K, Igusa invariants in K>
end function;

// Cheap screening variant: returns ONLY the absolute Igusa invariants over K (the
// period computation, without the Mestre descent). Use this to filter many x by
// good reduction at 5 before paying for a curve model.
function Genus2Invariants(x : K := Rationals(), Embedding := 0, Prec := 300,
                          Trials := 100, Verbose := false)
    CC := ComplexFieldExtra(Prec);
    tau, _, emb := x_to_tau_fast(x, CC : K := K, Embedding := Embedding,
                                 trials := Trials, verbose := Verbose);
    return IgusaInvariantsInK(tau, K, emb);
end function;
