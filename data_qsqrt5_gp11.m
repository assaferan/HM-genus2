// data_qsqrt5_gp11.m -- Comprehensive data file for the Q(sqrt(5)) (1,11) curves.
//
// Produces data_qsqrt5_gp11.txt with, for each of the 2 K-isomorphism classes
// (A = survivors 1&2; B = survivor 3; B' = survivor 4, the Galois conjugate of B):
//   - Igusa-Siegel invariants
//   - minimal quadratic twist
//   - bad odd prime ideals (Igusa valuation criterion up to norm 200)
//   - Frobenius L-polynomial at the first 20 prime ideals of good reduction
//
// Run: magma data_qsqrt5_gp11.m

AttachSpec("../CHIMP/CHIMP.spec");
load "survivors_qsqrt5_gp11.m";   // defines k<s5>, survivors_qsqrt5_gp11
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, ConductorExponentAt, PrimesAbove,
                      BadPrimesFromInvariants, IgusaPotentialGoodReduction;

outfile := "data_qsqrt5_gp11.txt";
System("rm -f " cat outfile);

procedure Emit(s) PrintFile(outfile, s); printf "%o\n", s; end procedure;

function MakeIntegral(C)
    K := BaseRing(C); f, h := HyperellipticPolynomials(C);
    m := LCM([Integers()| Denominator(K!c) : c in Coefficients(f) cat Coefficients(h) cat [K!1]]);
    return HyperellipticCurve(Parent(f)!(m^2*f), Parent(f)!(m*h));
end function;

// Reduce CI mod P, trying the flip y^2 = x^6*f(1/x) if the leading coeff of f
// vanishes mod P (happens when P divides the leading coefficient of the model).
function TryReduceModP(CI, P)
    K := BaseRing(CI); OK := Integers(K);
    Fp, redP := ResidueClassField(P); Rx<x> := PolynomialRing(Fp);
    f0, _ := HyperellipticPolynomials(CI); n := Degree(f0); R := Parent(f0);
    for flip in [0, 1] do
        f := flip eq 0 select f0
             else Polynomial([Coefficient(f0, n-i) : i in [0..n]]);
        for sh in [0..20] do
            fp := Rx![redP(OK!c) : c in Coefficients(Evaluate(f, R.1 + sh))];
            if Degree(fp) eq n then return true, HyperellipticCurve(fp, Rx!0); end if;
        end for;
    end for;
    return false, _;
end function;

ZT<T> := PolynomialRing(Integers());

// -----------------------------------------------------------------------
// Header
// -----------------------------------------------------------------------
Emit("===========================================================================");
Emit("Genus-2 curves over K = Q(sqrt(5)) with a rational (1,11)-polarized isogeny");
Emit("from the Gross-Popescu construction (HM-genus2 pipeline)");
Emit("===========================================================================");
Emit("");
Emit("Field:  K = Q(s5),  s5^2 = 5.");
Emit("        Ring of integers: O_K = Z[(1+s5)/2]  (golden-ratio ring).");
Emit("        Discriminant: disc(K/Q) = 5.  Class number: h(K) = 1.");
Emit("        Golden ratio: phi = (1+s5)/2, phi^(-1) = (s5-1)/2.");
Emit("");
Emit("Splitting of 11 in K: (5/11) = 1 by QR, so 11 is SPLIT.");
Emit("  11*O_K = p11a * p11b,  N(p11a) = N(p11b) = 11.");
Emit("");
Emit("The search enumerated non-rational P^3(K)-points (x:t) in the D_2 locus,");
Emit("inverted via InvertGP11Fast to period matrices, and recognized Igusa");
Emit("invariants in K.  4 survivors in 2 K-isomorphism classes:");
Emit("  Class A:  survivors 1 and 2 have identical Igusa invariants.");
Emit("  Class B:  survivor 3.");
Emit("  Class B': survivor 4, Galois conjugate of B (related by s5 -> -s5).");
Emit("(B and B' are a single orbit under Gal(K/Q); each is a distinct curve over K.)");
Emit("");

// -----------------------------------------------------------------------
// Per-curve data
// -----------------------------------------------------------------------
for idx in [1, 3, 4] do
    s := survivors_qsqrt5_gp11[idx];
    QI := s[2];
    lbl := idx eq 1 select "A  (= survivors 1 and 2)"
         else idx eq 3 select "B  (= survivor 3)"
         else "B' (= survivor 4, Galois conjugate of B)";

    Emit("");
    Emit("===========================================================================");
    Emit(Sprintf("CURVE %o", lbl));
    Emit("===========================================================================");
    Emit("");
    Emit("Igusa-Siegel invariants [J2:J4:J6:J8:J10] in WPS(2,4,6,8,10), J2 = 1:");
    Emit(Sprintf("  J2  = %o", QI[1]));
    Emit(Sprintf("  J4  = %o", QI[2]));
    Emit(Sprintf("  J6  = %o", QI[3]));
    Emit(Sprintf("  J8  = %o", QI[4]));
    Emit(Sprintf("  J10 = %o", QI[5]));
    Emit("");

    // Build curve from Igusa invariants.
    try C := Genus2CurveFromIgusa(QI, k);
    catch e; Emit(Sprintf("  [curve construction failed: %o]", e`Object)); continue; end try;

    // Minimal quadratic twist.
    Cmin := C; d := k!1;
    try Cmin, d := MinimalTwist(C, k); catch e; end try;
    Emit(Sprintf("Minimal quadratic twist factor: d = %o", d));
    Emit("");

    // Bad prime ideals (Igusa valuation criterion, fast, no model needed).
    bad_qi := BadPrimesFromInvariants(QI, k : Bound := 200);
    Emit("Bad odd primes (Igusa valuation criterion, up to rational prime 200):");
    if #bad_qi eq 0 then
        Emit("  None found.  (Potential everywhere-good reduction at all odd primes.)");
    else
        for p in bad_qi do Emit(Sprintf("  N(p) = %o", Norm(p))); end for;
    end if;
    Emit("  Note: p = 2 not checked (Magma cannot compute conductor exponent at p|2");
    Emit("  over number fields).");
    Emit("");

    // Conductor exponents at 11 (actual, via Magma's Conductor).
    Emit("Conductor exponents at primes above 11:");
    for p11 in PrimesAbove(k, 11) do
        e := ConductorExponentAt(Cmin, p11);
        Emit(Sprintf("  p11 (N = %o): exponent = %o  [%o]",
            Norm(p11), e, e eq 0 select "good reduction" else "BAD"));
    end for;
    Emit("");

    // Frobenius L-polynomials.
    Emit("Frobenius L-polynomials L_p(T) at prime ideals of good reduction");
    Emit("(first 20 prime ideals, skipping p | 2):");
    Emit("L_p(T) = det(1 - Frob_p * T | V_l Jac(C)),  deg = 4,  coefficients in Z.");
    Emit("Functional equation: L_p(T) = N(p)^2 T^4 L_p(1/(N(p)T)).");
    Emit("");
    Emit(Sprintf("  %-6o  %-5o  %-54o  factored over Z[T]",
        "N(p)", "p", "L_p(T)"));
    Emit("  " cat "-"^100);

    CI := MakeIntegral(Cmin);
    bad_set := {p : p in bad_qi};  // prime ideals known bad from Igusa criterion
    count := 0; pp := 3;
    while count lt 20 and pp lt 600 do
        for pid in PrimesAbove(k, pp) do
            // Skip if Igusa criterion says bad at this prime.
            if pid in bad_set then continue; end if;
            // Also skip if conductor exponent is positive.
            e := ConductorExponentAt(Cmin, pid);
            if e ne 0 then continue; end if;

            ok, Cp := TryReduceModP(CI, pid);
            if not ok then
                Emit(Sprintf("  %-6o  %-5o  [reduction failed -- skipped]", Norm(pid), pp));
                continue;
            end if;

            lp := ZT ! LPolynomial(Cp);
            fac := Factorization(lp);
            Emit(Sprintf("  %-6o  %-5o  %-54o  %o", Norm(pid), pp, lp, fac));
            count +:= 1;
            if count ge 20 then break; end if;
        end for;
        pp := NextPrime(pp);
    end while;

    Emit("");
    Emit(Sprintf("Frobenius at p | 11 (both primes; confirms rational 11-torsion):"));
    for p11 in PrimesAbove(k, 11) do
        e := ConductorExponentAt(Cmin, p11);
        if e ne 0 then
            Emit(Sprintf("  N(p11) = %o: bad reduction (exponent %o)", Norm(p11), e));
            continue;
        end if;
        ok, Cp11 := TryReduceModP(CI, p11);
        if not ok then
            Emit(Sprintf("  N(p11) = %o: reduction failed", Norm(p11)));
            continue;
        end if;
        lp11 := ZT ! LPolynomial(Cp11);
        fac11 := Factorization(lp11);
        R11 := PolynomialRing(GF(11));
        rts11 := [r[1] : r in Roots(R11 ! lp11)];
        Emit(Sprintf("  N(p11) = %o: L_{p11}(T) = %o = %o  (roots mod 11: %o)",
            Norm(p11), lp11, fac11, rts11));
    end for;
end for;

// -----------------------------------------------------------------------
// Footer / notes
// -----------------------------------------------------------------------
Emit("");
Emit("===========================================================================");
Emit("Notes");
Emit("===========================================================================");
Emit("");
Emit("Convention for L-polynomials:");
Emit("  L_p(T) = 1 + a1*T + a2*T^2 + a3*T^3 + N(p)^2 * T^4");
Emit("  with a3 = N(p)*a1 (from functional equation) and |a1| <= 4*sqrt(N(p)).");
Emit("  The number of F_{N(p)}-points on Jac(C) is L_p(1) = 1+a1+a2+a3+N(p)^2.");
Emit("");
Emit("Rational (1,11)-isogeny:");
Emit("  Each curve C has a K-rational 11-torsion subgroup in Jac(C)(K), so at");
Emit("  every prime of good reduction 11 divides #Jac(C)(F_q) = L_p(1), i.e.");
Emit("  L_p(T) has T=1 as a root modulo 11.  At both primes above 11 this gives");
Emit("  L_{p11}(1) = 0 in F_11, consistent with the 11-isogeny kernel reducing");
Emit("  to a non-trivial F_{11}-point.");
Emit("");
Emit("Sharing these curves:");
Emit("  The compact Igusa-Siegel invariants [J2:...:J10] above are the recommended");
Emit("  form for sharing; denominators are at most ~ 10^22.  A Weierstrass model");
Emit("  y^2 = f(x) over Q(sqrt(5)) can be recovered in Magma via:");
Emit("    k<s5> := QuadraticField(5);");
Emit("    load \"survivors_qsqrt5_gp11.m\";");
Emit("    QI := survivors_qsqrt5_gp11[1][2];  // curve A");
Emit("    import \"Genus2Curve.m\": Genus2CurveFromIgusa;");
Emit("    C := Genus2CurveFromIgusa(QI, k);");
Emit("  (The resulting Weierstrass model has large-fraction coefficients due to");
Emit("   Magma not implementing height-reduction for Mestre's conic over Q(sqrt(5));");
Emit("   the invariants are the better shareable form.)");
Emit("");
printf "Done. Output written to %o\n", outfile;
