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

// Local conductor exponent of C at p (an integer prime for Q, a prime ideal for K).
// Returns -1 if Magma cannot compute it (number-field fibre blowups at p | 2).
function ConductorExponentAt(C, p)
    try return Conductor(C, p); catch e; return -1; end try;
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
