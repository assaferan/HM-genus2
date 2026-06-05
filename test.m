// Test of the full pipeline x in P^3 -> genus-2 curve C.
//
// The HM construction produces a (1,5)-polarized abelian surface A_x; the
// reconstructed Jacobian Jac(C) is linked to it by a CYCLIC 5-ISOGENY. So Jac(C)
// carries a rational cyclic 5-isogeny -- a Galois-STABLE subgroup of order 5 --
// which is STRICTLY WEAKER than a rational 5-torsion point: Frobenius may act on
// that order-5 subgroup by any scalar in F_5^*, so its points need not be rational.
// (Hence the old `#TorsionSubgroup mod 5 eq 0` check tested the wrong thing.)
//
// Correct, computable necessary condition: a Galois-stable line in J[5] forces, at
// every prime p of good reduction (p != 5), Frobenius to have an eigenvalue in F_5
// on J[5] -- equivalently the mod-5 reduction of the L-polynomial of C/F_p has a
// root in F_5. We check this over many primes. It is a strong discriminator: a
// curve WITHOUT the 5-isogeny fails it at some small prime (a random genus-2 curve
// typically fails within the first few primes).
//
// NOTE: load Genus2Curve.m first so Genus2Curve is in scope:  magma Genus2Curve.m test.m

// Returns: <holds, #primes checked, first offending prime (0 if none)>.
function HasRational5IsogenyEvidence(C : Bound := 200, MinPrimes := 20)
    F5 := GF(5); R5 := PolynomialRing(F5);
    disc := Integers() ! (Numerator(Discriminant(C)) * Denominator(Discriminant(C)));
    checked := 0;
    for p in PrimesInInterval(7, Bound) do
        if p eq 5 or disc mod p eq 0 then continue; end if;   // skip p=5 and bad reduction
        Lp := R5 ! LPolynomial(ChangeRing(C, GF(p)));         // mod-5 L-polynomial of C/F_p
        if not HasRoot(Lp) then
            return false, checked, p;                          // no F_5-eigenvalue at p
        end if;
        checked +:= 1;
    end for;
    return checked ge MinPrimes, checked, 0;
end function;

procedure test()
    x := [1,2,3,4];
    C := Genus2Curve(x);
    printf "Reconstructed curve: %o\n", C;
    holds, checked, badp := HasRational5IsogenyEvidence(C);
    printf "Rational 5-isogeny evidence over %o good primes: %o\n", checked, holds;
    error if not holds,
        Sprintf("no rational 5-isogeny: Frobenius has no F_5-eigenvalue on J[5] at p=%o", badp);
    print "PASS: Jac(C) admits a rational cyclic 5-isogeny.";
end procedure;

test();
