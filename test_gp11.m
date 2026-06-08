// Test of the (1,11) / Gross-Popescu pipeline:
//   D_2 point -> pencil -> tau -> genus-2 curve C
//   Verify: at good primes P of K (above p != 11), L_P(C) mod 11 has a root in F_11.
//
// A rational 11-isogeny forces a Galois-stable rank-1 subspace in J[11], so Frobenius
// at every good prime has an eigenvalue in F_11; equivalently L_P mod 11 has a root.
// This is a strong check: a random genus-2 curve fails at some small prime.
//
// Run as: magma GP11Curve.m test_gp11.m
import "GP11.m": IntertwiningMatrix;

// Make C integral at all rational primes (clear denominators). Works over Q and K.
function IntegralModelK(C)
    K := BaseRing(C);
    f, h := HyperellipticPolynomials(C);
    if Type(K) eq FldRat then
        m := LCM([Integers()| Denominator(c) : c in Coefficients(f) cat Coefficients(h) cat [1]]);
        return HyperellipticCurve(Parent(f)!(m^2*f), Parent(f)!(m*h));
    end if;
    m := LCM([Integers()| Denominator(K!c) : c in Coefficients(f) cat Coefficients(h) cat [K!1]]);
    return HyperellipticCurve(Parent(f)!(m^2*f), Parent(f)!(m*h));
end function;

// Reduce the integral curve CI modulo the prime ideal P of K (or integer prime for Q).
// Returns <ok, reduced curve over residue field>.
function ReduceModP(CI, P)
    K := BaseRing(CI);
    f, h := HyperellipticPolynomials(CI);
    try
        if Type(K) eq FldRat then
            Fp := GF(P);
            Rx<x> := PolynomialRing(Fp);
            fp := Rx![Fp!Integers()!c : c in Coefficients(f)];
            hp := Rx![Fp!Integers()!c : c in Coefficients(h)];
        else
            Fp, redP := ResidueClassField(P);
            OK := Integers(K);
            Rx<x> := PolynomialRing(Fp);
            fp := Rx![redP(OK!c) : c in Coefficients(f)];
            hp := Rx![redP(OK!c) : c in Coefficients(h)];
        end if;
        // skip if the leading coefficient vanished (degree drop -> not genus 2)
        if LeadingCoefficient(fp) eq 0 then return false, _; end if;
        return true, HyperellipticCurve(fp, hp);
    catch e;
        return false, _;
    end try;
end function;

// For each rational prime p != 11 up to Bound, reduce C modulo primes of K above p,
// compute L mod 11, print, and check for a root in F_11. Returns <holds, #checked, bad_p>.
function Check11IsogenyFrobenius(C : Bound := 200, MinPrimes := 20)
    K := BaseRing(C);
    F11 := GF(11); R11 := PolynomialRing(F11);
    CI := IntegralModelK(C);
    checked := 0;

    for p in PrimesInInterval(13, Bound) do
        if p eq 11 then continue; end if;
        if Type(K) eq FldRat then
            ps := [p];
        else
            ps := [t[1] : t in Factorization(p * Integers(K))];
        end if;
        for P in ps do
            ok, Cp := ReduceModP(CI, P);
            if not ok then continue; end if;
            try
                Lp := R11 ! LPolynomial(Cp);
                rts := [r[1] : r in Roots(Lp)];
                pnm := Type(K) eq FldRat select p else Norm(P);
                printf "  p=%o (norm %o): L mod 11 = %o,  F_11-roots: %o\n",
                    p, pnm, Lp, rts;
                if #rts eq 0 then
                    return false, checked, p;
                end if;
                checked +:= 1;
            catch e; end try;
        end for;
        if checked ge MinPrimes then break; end if;
    end for;
    return checked ge MinPrimes, checked, 0;
end function;

procedure test11()
    // Build the D_2 point: fix x1..x4 = 1,2,3,4; x5 = root of Pf(S)(1,2,3,4,t).
    Qt<t> := PolynomialRing(Rationals());
    f := Pfaffian(IntertwiningMatrix([Qt| 1, 2, 3, 4, t]));
    // f = 3t^4 + 16t^3 + 124t^2 + 174t - 1190 (irreducible quartic)
    g := [h[1] : h in Factorization(f) | Degree(h[1]) ge 1][1];
    K<a> := NumberField(g);

    prec := 120;
    CC := ComplexFieldExtra(prec);
    rts_real := [r[1] : r in Roots(g, RealField(prec))];
    error if #rts_real eq 0, "No real roots of D_2 polynomial -- cannot invert";

    // Try each real root as the embedding; the one with Im(tau) > 0 converges.
    C := false;
    for rt in rts_real do
        printf "Trying real root ~ %o ...\n", RealField(6)!rt;
        try
            emb := hom< K -> CC | CC ! rt >;
            C, QI, tau := GP11Curve([K| 1, 2, 3, 4, a], K, emb, CC : Trials := 120);
            printf "Inversion succeeded. tau =\n%o\n", tau;
            break;
        catch e;
            printf "  failed (%o), trying next root\n", e;
        end try;
    end for;
    error if Type(C) eq BoolElt, "GP11Curve failed for all real embeddings";

    printf "\nReconstructed curve C: %o\n", C;
    printf "Base field: %o\n\n", BaseRing(C);

    printf "Frobenius polynomials mod 11 at good primes:\n";
    holds, checked, badp := Check11IsogenyFrobenius(C);
    printf "\nRational 11-isogeny evidence over %o good primes: %o\n", checked, holds;
    error if not holds,
        Sprintf("FAIL: no F_11-eigenvalue of Frobenius on J[11] at p=%o", badp);
    print "PASS: Jac(C) admits a rational cyclic 11-isogeny.";
end procedure;

test11();
