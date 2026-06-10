// Reduction / conductor utilities for the reconstructed genus-2 curve C, over Q or a
// number field K. Two layers:
//   (a) a cheap, model-free filter from the Igusa invariants (potential good reduction
//       at p via the valuation criterion) -- usable directly on the invariants the
//       pipeline already produces, before/without building the curve model;
//   (b) the actual local conductor exponent via Magma's Conductor(C, p), which over a
//       number field works at all ODD primes (in particular p | 5) but is not yet
//       implemented at p | 2.

// Primes of O_K above the rational prime ell (just [ell] when K = Q).
function PrimesAbove(K, ell)
    if Type(K) eq FldRat then return [ell]; end if;
    return [ t[1] : t in Factorization(ell*Integers(K)) ];
end function;

// Igusa valuation criterion for POTENTIAL good reduction at p (residue char != 2):
// scaling J_{2i} -> u^{2i} J_{2i} can make all invariants integral with J10 a unit iff
//   v_p(J10)/10 <= v_p(J_{2i})/(2i)  for i = 1..4.
// (Scaling-invariant, so any weighted-projective representative of the invariants works,
// e.g. the normalized invariants the reconstruction returns.) Jinv = [J2,J4,J6,J8,J10].
function IgusaPotentialGoodReduction(Jinv, p)
    wts := [2,4,6,8,10];
    base := Valuation(Jinv[5], p) / 10;
    for i in [1..4] do
        if Jinv[i] ne 0 and Valuation(Jinv[i], p)/wts[i] lt base then return false; end if;
    end for;
    return true;
end function;

// Potential good reduction at every prime above 5 (model-free, from the invariants).
function PotentialGoodReductionAt5(Jinv, K)
    return &and[ IgusaPotentialGoodReduction(Jinv, p) : p in PrimesAbove(K, 5) ];
end function;

// All prime ideals in the support of the Igusa invariants (away from 2).
// These are the ONLY primes that can fail the valuation criterion: any prime p
// with v_p(J_i) = 0 for all i automatically satisfies v_p(J_{2i})/(2i) >= v_p(J10)/10.
function JinvSupport(Jinv, K)
    OK := Integers(K);
    prime_set := {};
    for j in Jinv do
        if j eq 0 then continue; end if;
        if Type(K) eq FldRat then
            n := Integers() ! Numerator(j); d := Integers() ! Denominator(j);
            for fac in Factorization(Abs(n) * d) do
                if fac[1] ne 2 then Include(~prime_set, fac[1]); end if;
            end for;
        else
            for fac in Factorization(j * OK) do
                p := fac[1];
                if Norm(p) mod 2 ne 0 then Include(~prime_set, p); end if;
            end for;
        end if;
    end for;
    return prime_set;
end function;

// Bad prime ideals from the Igusa valuation criterion.
// Exact (no bound needed): only primes in the support of the invariants can be bad.
function BadPrimesFromInvariants(Jinv, K)
    bad := [];
    for p in JinvSupport(Jinv, K) do
        if not IgusaPotentialGoodReduction(Jinv, p) then Append(~bad, p); end if;
    end for;
    return bad;
end function;

// A "conductor size" proxy for ranking survivors: the norm of the product of the
// bad prime ideals (from the Igusa criterion). Smaller is nicer.
// Returns <size, bad primes>.
function ConductorSizeProxy(Jinv, K)
    bad := BadPrimesFromInvariants(Jinv, K);
    return &*([1] cat [ Type(K) eq FldRat select p else Norm(p) : p in bad ]), bad;
end function;

// Local conductor exponent of C at p (an integer prime for Q, a prime ideal for K).
// Returns -1 if Magma cannot compute it (number-field fibre blowups at p | 2).
function ConductorExponentAt(C, p)
    try return Conductor(C, p); catch e; return -1; end try;
end function;

// Integral model of C: clear denominators of the hyperelliptic polynomials. y^2+h y=f
// is isomorphic to Y^2 + (m h) Y = m^2 f (Y = m y), with m the common denominator -- so
// the reduction (hence Euler factors/conductor) is unchanged but the model is integral.
function IntegralModelOf(C)
    f, h := HyperellipticPolynomials(C);
    m := LCM([ Integers() | Denominator(c) : c in Coefficients(f) cat Coefficients(h) cat [1] ]);
    return HyperellipticCurve(Parent(f)!(m^2*f), Parent(f)!(m*h));
end function;

// Frobenius polynomial (Euler factor) of C at a good prime p, as an integer polynomial.
// Uses an integral model (EulerFactor needs a p-integral model; the reconstructed curve
// usually is not). Degenerates if C has bad reduction at p.
function FrobeniusPolynomial(C, p)
    return PolynomialRing(Integers()) ! EulerFactor(IntegralModelOf(C), p);
end function;

// Minimal quadratic twist (best effort): twist C by the product of generators of the
// odd prime ideals (norm via rational primes <= Bound) where C has additive bad
// reduction -- those are primes of potential good reduction (invariant criterion), and
// twisting by a uniformizer there toggles them to good without disturbing other primes
// (the generator is a unit elsewhere; needs class number 1). The result is good at all
// such primes when the potential-good reduction is of quadratic type (verify via the
// conductor). Returns <twisted curve, twist element>.
function MinimalTwist(C, K : Bound := 50)
    if Type(K) eq FldRat then
        bad := [ ell : ell in PrimesUpTo(Bound) | ell ne 2 and ConductorExponentAt(C, ell) gt 0 ];
        d := &*([Integers()| 1] cat bad);
        return (d eq 1 select C else QuadraticTwist(C, d)), d;
    end if;
    gens := [K| ];
    for ell in PrimesUpTo(Bound) do
        if ell eq 2 then continue; end if;
        for p in PrimesAbove(K, ell) do
            if ConductorExponentAt(C, p) gt 0 then
                ok, pi := IsPrincipal(p);
                error if not ok, "MinimalTwist: non-principal prime (class number > 1)";
                Append(~gens, K!pi);
            end if;
        end for;
    end for;
    d := &*([K| 1] cat gens);
    return (d eq 1 select C else QuadraticTwist(C, d)), d;
end function;

// True iff C has (actual) good reduction at every prime above 5 -- conductor exponent
// 0 there. Definitive (p | 5 is odd, so Conductor(C,p) is computable over K).
function GoodReductionAt5(C)
    return &and[ ConductorExponentAt(C, p) eq 0 : p in PrimesAbove(BaseRing(C), 5) ];
end function;

// Print a reduction summary: status at 5, and the conductor (exponents per bad prime;
// the 2-part is flagged when Magma cannot compute it over a number field).
procedure ReductionReport(C)
    K := BaseRing(C);
    J := IgusaInvariants(C);
    printf "Reduction at primes above 5:\n";
    for p in PrimesAbove(K, 5) do
        ce := ConductorExponentAt(C, p);
        printf "  p (norm %o): potential-good=%o  exponent=%o  -> %o\n",
            Type(K) eq FldRat select p else Norm(p),
            IgusaPotentialGoodReduction(J, p), ce,
            ce eq 0 select "GOOD" else (ce eq -1 select "not computable" else "bad");
    end for;
    if Type(K) eq FldRat then
        printf "Conductor = %o\n", Factorization(Conductor(C));
    else
        bad := [ t[1] : t in Factorization(Discriminant(C)*Integers(K)) ];
        printf "Conductor exponents (bad primes; norm:exponent, -1 = not computable at 2):\n";
        for p in bad do printf "  norm %o : %o\n", Norm(p), ConductorExponentAt(C, p); end for;
    end if;
end procedure;
