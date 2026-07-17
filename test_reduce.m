// Diagnose TryReduceModP failures: print original (non-integral) P-valuations.
AttachSpec("../CHIMP/CHIMP.spec");
load "survivors.m";
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, PrimesAbove;

k<s2> := QuadraticField(2);
OK := Integers(k);
ZT<T> := PolynomialRing(Integers());

// Reduce a K-rational (possibly non-integral) element mod P:
// For K = Q(sqrt(2)) and P inert over ell, elements are a + b*s2 with a,b ∈ Q.
// P-integral iff v_ell(a) >= 0 and v_ell(b) >= 0 (for the ell-part of a,b).
function ReduceKElt(c, redP, Fp, ell)
    coords := Eltseq(k!c);  // [a_rat, b_rat] with c = a + b*s2
    function RedQ(q)
        n := Integers()!Numerator(q); d := Integers()!Denominator(q);
        error if d mod ell eq 0, "element not P-integral";
        return Fp!n * (Fp!d)^(-1);
    end function;
    return RedQ(coords[1]) + RedQ(coords[2]) * redP(OK.1);
end function;

// New TryReduceModP: takes C (possibly non-integral) and handles K-rational coefficients.
function TryReduceModP(C, P)
    K := BaseRing(C); OK2 := Integers(K);
    Fp, redP := ResidueClassField(P); Rx<t> := PolynomialRing(Fp);
    ell := Characteristic(Fp);
    f0, _ := HyperellipticPolynomials(C); n := Degree(f0); Rk<x> := Parent(f0);
    ww := K.1;

    // Residue classes in O_K for shifts
    res := [K| a + b*ww : a in [0..ell-1], b in [0..ell-1]];

    // P-adic valuation of a K-element (handles non-integral, for inert ell)
    function VP(c)
        if K!c eq 0 then return Infinity(); end if;
        coords := Eltseq(K!c);
        va := Valuation(Numerator(coords[1]), ell) - Valuation(Denominator(coords[1]), ell);
        vb := Valuation(Numerator(coords[2]), ell) - Valuation(Denominator(coords[2]), ell);
        return Min(va, vb);
    end function;

    // Scale polynomial to P-minimal form: find e so that
    // all a_i * ell^{(i-n)*e} are P-integral and min valuation = 0.
    // Equivalently: e = floor(min_{i<n} v_P(a_i)/(n-i)).
    // After scaling: new a_i = a_i * ell^{(i-n)*e}.
    function ScaleMinimal(f)
        nn := Degree(f);
        min_rat := Infinity();
        for i in [0..nn-1] do
            c := Coefficient(f, i);
            if K!c ne 0 then
                v := VP(c);
                r := Rationals()!v / (nn - i);
                if r lt min_rat then min_rat := r; end if;
            end if;
        end for;
        if min_rat eq Infinity() then return f; end if;
        e := Floor(min_rat);
        new_coeffs := [K!Coefficient(f,i) * ell^((i-nn)*e) : i in [0..nn]];
        return Polynomial(new_coeffs);
    end function;

    // Try shifting polynomial poly (with K-rational, P-integral coefficients), reduce mod P.
    function TryShifts(poly)
        d := Degree(poly);
        for sh in res do
            g := Evaluate(poly, Rk.1 + sh);
            cffs := [K!Coefficient(g, i) : i in [0..d]];
            if exists{c : c in cffs | VP(c) lt 0} then continue; end if;
            ok2 := true; fp_coeffs := [];
            for c in cffs do
                try Append(~fp_coeffs, ReduceKElt(c, redP, Fp, ell));
                catch e; ok2 := false; break; end try;
            end for;
            if not ok2 then continue; end if;
            fp := Polynomial(fp_coeffs);
            if Degree(fp) eq d and Discriminant(fp) ne 0 then
                return true, HyperellipticCurve(fp, Rx!0);
            end if;
        end for;
        return false, _;
    end function;

    for use_flip in [false, true] do
        f_base := use_flip
            select Polynomial([K!Coefficient(f0, n-i) : i in [0..n]])
            else f0;
        f := ScaleMinimal(f_base);
        ok2, Cp := TryShifts(f); if ok2 then return ok2, Cp; end if;
        // Mobius: x -> alpha + 1/X, leading coeff = f_scaled(alpha)
        for alpha in res do
            fa := Evaluate(f, alpha);
            if VP(fa) gt 0 then continue; end if;
            g := &+ [K!Coefficient(f, i) * x^(n-i) * (alpha*x + 1)^i : i in [0..n]];
            ok2, Cp := TryShifts(g); if ok2 then return ok2, Cp; end if;
        end for;
    end for;
    return false, _;
end function;

for idx in [2, 12, 13] do
    s := survivors[idx];
    QI := s[2];
    printf "=== Class %o: x=%o ===\n", idx, s[1];
    C := Genus2CurveFromIgusa(QI, k);
    Cmin, d := MinimalTwist(C, k);
    f_raw, _ := HyperellipticPolynomials(Cmin);  // raw, possibly non-integral

    // Print P5-valuations of the RAW (non-integral) polynomial
    p5 := PrimesAbove(k, 5)[1];
    p3 := PrimesAbove(k, 3)[1];
    printf "  P5-valuations of Cmin polynomial: ";
    for i in [0..Degree(f_raw)] do
        c := Coefficient(f_raw, i);
        if k!c eq 0 then printf "inf "; continue; end if;
        coords := Eltseq(k!c);
        va := Valuation(Numerator(coords[1]), 5) - Valuation(Denominator(coords[1]), 5);
        vb := Valuation(Numerator(coords[2]), 5) - Valuation(Denominator(coords[2]), 5);
        printf "%o ", Min(va,vb);
    end for; printf "\n";

    printf "  Testing p5 (N=25, using Cmin directly): ";
    ok2, Cp := TryReduceModP(Cmin, p5);
    if ok2 then
        lp := ZT ! LPolynomial(Cp); printf "SUCCESS! L = %o\n", lp;
    else printf "FAIL\n"; end if;

    printf "  Testing p3 (N=9, using Cmin directly): ";
    ok2, Cp := TryReduceModP(Cmin, p3);
    if ok2 then
        lp := ZT ! LPolynomial(Cp); printf "SUCCESS! L = %o\n", lp;
    else printf "FAIL\n"; end if;
end for;
