// Genus-2 curve reconstruction for a point x in P^3(K) over a real quadratic field
// K = Q(sqrt(d)) -- the number-field analogue of Genus2Curve.m (the Q case).
//
// Pipeline (same shape as the Q case, with two changes):
//   1. the HM surface A_x is built over K, so its defining quintics live over K;
//   2. a FIXED real embedding sqrt(d) -> EmbSign*|sqrt(d)| in CC turns those quintics
//      into complex polynomials -- the numeric inversion is then identical;
//   3. the curve's absolute Igusa invariants are recognized as elements of K (degree
//      <= 2) under the SAME embedding, and the curve is reconstructed and DESCENDED
//      to K via Mestre's algorithm (HyperellipticCurveFromIgusaClebsch).
//
// The embedding choice EmbSign = +1 / -1 selects between the Galois-conjugate inputs
// (sqrt(d) -> +1.414... vs -1.414...) and hence conjugate curves; use it consistently.
//
// Worked example (x = [1, sqrt2, 3, 4] over Q(sqrt2)), absolute Igusa invariants:
//   I1 = 1
//   I2 = (342163089294935*sqrt2 - 483535770283008)/14874181170436
//   ... (degree-2 over Q(sqrt2); see README). The curve descends to Q(sqrt2).
AttachSpec("../CHIMP/CHIMP.spec");
import "InversionFast.m": FindTauFast;
import "Inversion.m": verify_tau;

// Recognize a complex value as an element of K (degree <= [K:Q]) under embedding emb.
// Returns <element of K, recognition residual>.
function RecognizeInField(z, K, emb)
    mp := MinimalPolynomial(z, Degree(K));
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
// Same theta-derivative computation as Genus2CurveFromTau in Genus2Curve.m.
function IgusaNumericFromTau(tau)
    CC := BaseRing(tau);
    P := HorizontalJoin(tau, IdentityMatrix(CC, 2));
    taunew, gamma := ReduceSmallPeriodMatrix(tau);
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

// Reconstruct the curve over K from a period matrix tau, recognizing the Igusa
// invariants in K (embedding emb) and descending via Mestre. Best-effort reduction.
function CurveOverFieldFromTau(tau, K, emb : Reduce := true)
    ICCn := IgusaNumericFromTau(tau);
    QI := [K| ];
    for z in ICCn do
        el, res := RecognizeInField(z, K, emb);
        error if res gt 1e-20,
            "CurveOverFieldFromTau: Igusa invariant not recognized in K; increase Prec";
        Append(~QI, el);
    end for;
    // HyperellipticCurveFromIgusaInvariants may return a model over an extension of K;
    // its Igusa-Clebsch invariants still lie in K, and Mestre's algorithm descends to K
    // whenever the Mestre conic has a K-rational point.
    Cext := HyperellipticCurveFromIgusaInvariants(QI);
    IC := [K| c : c in IgusaClebschInvariants(Cext) ];
    C := HyperellipticCurveFromIgusaClebsch(IC);
    if Reduce then
        try C := ReducedMinimalWeierstrassModel(C); catch e
            try C := ReducedModel(C); catch e2 ; end try;
        end try;
    end if;
    return C, QI;
end function;

// Full pipeline: x in P^3(K) -> genus-2 curve over K (when it descends).
//   K            -- a real quadratic field Q(sqrt(d)) (e.g. K<s2> := QuadraticField(2))
//   x            -- sequence of 4 elements of K
//   Prec         -- working precision for the polish / invariant recognition
//   SearchPrec   -- reduced precision for the random-restart basin search
//   Trials       -- max restarts (only ~10% land in the true basin)
//   EmbSign      -- +1 or -1: which real embedding of sqrt(d) to use
function Genus2CurveQuadratic(x, K : Prec := 300, SearchPrec := 40, Trials := 200,
                              EmbSign := 1, BasinThresh := 1e-6, Verbose := false)
    d := Rationals() ! (K.1^2);                  // K = Q(sqrt(d))
    // --- HM surface over K and its defining quintics ---
    P := ProjectiveSpace(K, 4);
    F := HorrocksMumfordBundle(P);
    M := GlobalSectionSubmodule(F);
    secs := [M.i : i in [1..Rank(M)] | Grading(M)[i] eq 0];
    assert #secs eq 4;
    s := &+[x[i]*secs[i] : i in [1..4]];
    A := Scheme(Ambient(ZeroSubscheme(F, s)), Saturation(Ideal(ZeroSubscheme(F, s))));
    qs := [b : b in MinimalBasis(Ideal(A)) | Degree(b) eq 5];
    assert #qs eq 3;
    P5k := Parent(qs[1]);

    // coerce the K-quintics to CC via the chosen real embedding of sqrt(d)
    coerce := function(CC)
        emb := hom< K -> CC | EmbSign*Sqrt(CC!d) >;
        P5CC := PolynomialRing(CC, 5);
        cm := hom< P5k -> P5CC | emb, [P5CC.j : j in [1..5]] >;
        qsCC := [ cm(q) : q in qs ];
        return qsCC, [ [ Derivative(q, l) : l in [1..5] ] : q in qsCC ], emb;
    end function;

    // --- Phase A: locate tau at reduced precision ---
    CCs := ComplexFieldExtra(SearchPrec); iis := CCs.1;
    qsS, dqsS := coerce(CCs);
    zs := [ [CCs| (Random(-50,50)+iis*Random(10,40))/100,
                  (Random(-50,50)+iis*Random(10,40))/100 ] : kk in [1..6] ];
    best := [];
    for trial in [1..Trials] do
        tau0 := [ (Random(-200,200)+iis*Random(10,120))/100,
                  (Random(-150,150)+iis*Random(-50,50))/100,
                  (Random(-250,250)+iis*Random(10,250))/100 ];
        ts, rs := FindTauFast(qsS, dqsS, tau0, zs, CCs : iters := 50, verbose := Verbose);
        pd := Imaginary(ts[1]) gt 0 and Imaginary(ts[1])*Imaginary(ts[3]) - Imaginary(ts[2])^2 gt 0;
        if pd and rs lt BasinThresh then best := ts; break; end if;
    end for;
    error if #best eq 0, "Genus2CurveQuadratic: no tau found; increase Trials";

    // --- Phase B: polish at full precision, then reconstruct over K ---
    CC := ComplexFieldExtra(Prec);
    qsH, dqsH, emb := coerce(CC);
    tau0H := [ CC!Real(t) + CC.1*(CC!Imaginary(t)) : t in best ];
    zsH := [ [CC| (Random(-50,50)+CC.1*Random(10,40))/100,
                  (Random(-50,50)+CC.1*Random(10,40))/100 ] : kk in [1..6] ];
    tauv := FindTauFast(qsH, dqsH, tau0H, zsH, CC : iters := 40, coarse := false,
                        lm_floor := 0, tol_pow := Prec-10, verbose := Verbose);
    verify_tau(qsH, tauv, CC);
    return CurveOverFieldFromTau(SymmetricMatrix(tauv), K, emb);
end function;
