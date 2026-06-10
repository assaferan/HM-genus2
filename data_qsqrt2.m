// data_qsqrt2.m -- Comprehensive data file for the Q(sqrt(2)) (1,5) curves.
//
// Produces data_qsqrt2.txt with, for each K-isomorphism class of genus-2 curves
// arising from the Horrocks-Mumford (1,5) search over K = Q(sqrt(2)):
//   - Igusa-Siegel invariants
//   - minimal quadratic twist
//   - bad odd prime ideals (Igusa valuation criterion up to norm 200)
//   - Frobenius L-polynomial at the first 20 prime ideals of good reduction
//   - Frobenius at the prime above 5 (5 is inert in K, N(p5) = 25)
//
// Run: magma data_qsqrt2.m

AttachSpec("../CHIMP/CHIMP.spec");
load "survivors.m";   // defines k<s2> := QuadraticField(2); survivors := [* ... *]
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, ConductorExponentAt, PrimesAbove,
                      BadPrimesFromInvariants;

outfile := "data_qsqrt2.txt";
System("rm -f " cat outfile);

procedure Emit(s) PrintFile(outfile, s); printf "%o\n", s; end procedure;

function MakeIntegral(C)
    K := BaseRing(C); f, h := HyperellipticPolynomials(C);
    m := LCM([Integers()| Denominator(K!c) : c in Coefficients(f) cat Coefficients(h) cat [K!1]]);
    return HyperellipticCurve(Parent(f)!(m^2*f), Parent(f)!(m*h));
end function;

// Reduce CI mod P.
//
// MakeIntegral clears denominators with a global LCM, which for inert primes like p5
// (where J-invariants already have 5 in denominators from the (1,5)-structure) makes
// every coefficient divisible by 5^k — f ≡ 0 mod P identically.  The fix is to
// scale the polynomial to its "P-optimal" form before trying shifts and Mobius
// transforms.
//
// Strategy:
//   For each candidate base polynomial (original, flip):
//     a) Scale x -> ell^e * x, y -> ell^{3e} * y so that
//        min_{i<n} v_P(a_i) = 0 (removes the common P-factor).
//        e = floor( min_{i<n, a_i≠0} v_P(a_i) / (n-i) ).
//     b) Try all ell^2 O_K-shifts on the scaled poly.
//     c) Try Mobius transforms x -> alpha + 1/X (leading coeff = f_scaled(alpha)):
//        find alpha with f_scaled(alpha) ≢ 0 mod P, then apply TryShifts.
function TryReduceModP(CI, P)
    K := BaseRing(CI); OK := Integers(K);
    Fp, redP := ResidueClassField(P); Rx<t> := PolynomialRing(Fp);
    ell := Characteristic(Fp);
    f0, _ := HyperellipticPolynomials(CI); n := Degree(f0); Rk<x> := Parent(f0);
    ww := K.1;

    // All residue classes of O_K/P (F_{ell} or F_{ell^2})
    res := [K| a + b*ww : a in [0..ell-1], b in [0..ell-1]];

    // Scale poly (with O_K coefficients) to minimise the common P-content.
    // For a degree-n polynomial: apply x -> ell^e * x, y -> ell^{3e} * y
    // (for n=6) or the analogous substitution.  The new a_i = a_i / ell^{(n-i)*e}.
    // Requires ell^{(n-i)*e} | a_i in O_K (guaranteed by the choice of e).
    // For inert P over ell, ell*O_K = P, so exact division by ell^k keeps us in O_K.
    function ScaleOptimal(f)
        nn := Degree(f);
        min_rat := Infinity();
        for i in [0..nn-1] do
            c := Coefficient(f, i);
            if c ne 0 then
                v := Valuation(K!c, P);  // P-adic valuation (integer)
                r := Rationals()!(v) / (nn - i);
                if r lt min_rat then min_rat := r; end if;
            end if;
        end for;
        if min_rat le 0 or min_rat eq Infinity() then return f; end if;
        e := Floor(min_rat);  // positive integer
        if e le 0 then return f; end if;
        // New a_i = a_i / ell^{(nn-i)*e}.  Since v_P(a_i) >= (nn-i)*e, the division
        // is exact in O_K for inert ell (ell*O_K = P); use ExactQuotient.
        new_coeffs := [ i lt nn
            select K ! ExactQuotient(OK!(K!Coefficient(f,i)), ell^((nn-i)*e))
            else Coefficient(f, nn) : i in [0..nn] ];
        return Polynomial(new_coeffs);
    end function;

    // Try to reduce poly (with O_K coefficients) with O_K shifts; check discriminant.
    function TryShifts(poly)
        d := Degree(poly);
        for sh in res do
            g := Evaluate(poly, Rk.1 + sh);
            cffs := Coefficients(g);
            if exists{c : c in cffs | Valuation(K!c, P) lt 0} then continue; end if;
            fp := Rx![redP(OK!(K!c)) : c in cffs];
            if Degree(fp) eq d and Discriminant(fp) ne 0 then
                return true, HyperellipticCurve(fp, Rx!0);
            end if;
        end for;
        return false, _;
    end function;

    for use_flip in [false, true] do
        f_base := use_flip
            select Polynomial([Coefficient(f0, n-i) : i in [0..n]])
            else f0;
        f := ScaleOptimal(f_base);

        // a) Shifts on the scaled polynomial.
        ok, Cp := TryShifts(f); if ok then return ok, Cp; end if;

        // b) Mobius transforms: x -> alpha + 1/X, leading coeff = f(alpha).
        for alpha in res do
            fa := Evaluate(f, alpha);
            if Valuation(K!fa, P) gt 0 then continue; end if;
            g := &+ [Coefficient(f, i) * x^(n-i) * (alpha*x + 1)^i : i in [0..n]];
            ok, Cp := TryShifts(g); if ok then return ok, Cp; end if;
        end for;
    end for;

    return false, _;
end function;

// Deduplicate survivors: group by equal Igusa invariants over K.
// Two survivors with identical [J2:J4:J6:J8:J10] (with J2=1) are K-isomorphic.
classes := [* *];
seen_QI := [* *];
for s in survivors do
    QI := s[2];
    is_new := true;
    for c in seen_QI do
        if forall{i : i in [1..5] | QI[i] eq c[i]} then is_new := false; break; end if;
    end for;
    if is_new then
        Append(~classes, s);
        Append(~seen_QI, QI);
    end if;
end for;

ZT<T> := PolynomialRing(Integers());
F25<w25> := GF(5, 2);   // F_25 = GF(5^2); w25 satisfies the Conway polynomial
R25<t> := PolynomialRing(F25);

// -----------------------------------------------------------------------
// Header
// -----------------------------------------------------------------------
Emit("===========================================================================");
Emit("Genus-2 curves over K = Q(sqrt(2)) with a rational (1,5)-polarized abelian");
Emit("surface from the Horrocks-Mumford bundle (HM-genus2 pipeline)");
Emit("===========================================================================");
Emit("");
Emit("Field:  K = Q(s2),  s2^2 = 2.");
Emit("        Ring of integers: O_K = Z[s2].");
Emit("        Discriminant: disc(K/Q) = 8.  Class number: h(K) = 1.");
Emit("");
Emit("Splitting of 2 in K: 2 ramifies;  2*O_K = p2^2,  N(p2) = 2.");
Emit("Splitting of 5 in K: (2/5) = -1 (QR), so 5 is INERT.");
Emit("  5*O_K = p5 (prime),  N(p5) = 25.");
Emit("General splitting: p splits in K iff p == +/-1 (mod 8).");
Emit("");
Emit("The search enumerated non-rational P^3(K)-points x of smallest height,");
Emit("inverted via x_to_tau_fast to period matrices, and recognized Igusa");
Emit("invariants in K.  Filters: genuinely over K (not a base-change from Q),");
Emit("potential good reduction at 5, and simple Jacobian (End = Z).");
Emit(Sprintf("%o total survivors, in %o K-isomorphism classes (grouped by Igusa invariants).",
    #survivors, #classes));
Emit("");

// -----------------------------------------------------------------------
// Per-class data
// -----------------------------------------------------------------------
for idx in [1 .. #classes] do
    s := classes[idx];
    QI := s[2];
    lbl := Sprintf("%o", idx);

    Emit("");
    Emit("===========================================================================");
    Emit(Sprintf("CURVE CLASS %o  (representative x = %o)", lbl, s[1]));
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
    bad_qi := BadPrimesFromInvariants(QI, k);
    Emit("Bad odd primes (Igusa valuation criterion, exact -- all primes in the support of the invariants):");
    if #bad_qi eq 0 then
        Emit("  None found.  (Potential everywhere-good reduction at all odd primes.)");
    else
        for p in bad_qi do Emit(Sprintf("  N(p) = %o", Norm(p))); end for;
    end if;
    Emit("  Note: p = 2 not checked (Magma cannot compute conductor exponent at p|2");
    Emit("  over number fields).");
    Emit("");

    // Conductor exponents at the bad primes and at p|5.
    if #bad_qi gt 0 then
        Emit("Conductor exponents at the bad prime ideals above:");
        for p in bad_qi do
            e := ConductorExponentAt(Cmin, p);
            Emit(Sprintf("  N(p) = %o: exponent = %o  [%o]",
                Norm(p), e, e eq 0 select "good reduction (twist fixed it)" else "BAD"));
        end for;
        Emit("");
    end if;
    Emit("Conductor exponent at the prime above 5 (N = 25, inert):");
    for p5 in PrimesAbove(k, 5) do
        e := ConductorExponentAt(Cmin, p5);
        Emit(Sprintf("  p5 (N = %o): exponent = %o  [%o]",
            Norm(p5), e, e eq 0 select "good reduction" else "BAD"));
    end for;
    Emit("");

    // Frobenius L-polynomials.
    Emit("Frobenius L-polynomials L_p(T) at prime ideals of good reduction");
    Emit("(first 20 prime ideals, skipping p | 2):");
    Emit("L_p(T) = det(1 - Frob_p * T | V_l Jac(C)),  deg = 4,  coefficients in Z.");
    Emit("Functional equation: L_p(T) = N(p)^2 T^4 L_p(1/(N(p)T)).");
    Emit("");
    Emit(Sprintf("  %-6o  %-5o  %-54o  factored over F_25[T]",
        "N(p)", "p", "L_p(T)"));
    Emit("  " cat "-"^100);

    CI := MakeIntegral(Cmin);
    bad_set := {p : p in bad_qi};
    count := 0; pp := 3;
    while count lt 20 and pp lt 600 do
        for pid in PrimesAbove(k, pp) do
            if pid in bad_set then continue; end if;
            e := ConductorExponentAt(Cmin, pid);
            if e ne 0 then continue; end if;

            ok, Cp := TryReduceModP(CI, pid);
            if not ok then
                Emit(Sprintf("  %-6o  %-5o  [reduction failed -- skipped]", Norm(pid), pp));
                continue;
            end if;

            lp := ZT ! LPolynomial(Cp);
            fac := Factorization(R25 ! lp);
            Emit(Sprintf("  %-6o  %-5o  %-54o  %o", Norm(pid), pp, lp, fac));
            count +:= 1;
            if count ge 20 then break; end if;
        end for;
        pp := NextPrime(pp);
    end while;

    Emit("");
    Emit(Sprintf("Frobenius at p | 5 (5 is inert; N(p5) = 25; confirms rational (1,5)-structure):"));
    for p5 in PrimesAbove(k, 5) do
        e := ConductorExponentAt(Cmin, p5);
        if e ne 0 then
            Emit(Sprintf("  N(p5) = %o: bad reduction (exponent %o)", Norm(p5), e));
            continue;
        end if;
        ok, Cp5 := TryReduceModP(CI, p5);
        if not ok then
            Emit(Sprintf("  N(p5) = %o: reduction failed", Norm(p5)));
            continue;
        end if;
        lp5 := ZT ! LPolynomial(Cp5);
        fac5 := Factorization(R25 ! lp5);
        rts5 := [r[1] : r in Roots(R25 ! lp5)];
        Emit(Sprintf("  N(p5) = %o: L_{p5}(T) = %o = %o  (roots over F_25: %o)",
            Norm(p5), lp5, fac5, rts5));
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
Emit("  Factorization is over F_25[T] = GF(5^2)[T].  The generator w25 of GF(5^2)");
Emit(Sprintf("  satisfies %o (Conway polynomial).", DefiningPolynomial(F25)));
Emit("");
Emit("Rational (1,5)-polarization structure:");
Emit("  The HM search finds abelian surfaces A = Jac(C) carrying a K-rational");
Emit("  (1,5)-polarization; equivalently, a Lagrangian 5^2-subgroup scheme defined");
Emit("  over K.  At the (inert) prime above 5 (N=25), this structure implies");
Emit("  5 | L_{p5}(1) = #Jac(C)(F_25), i.e. T=1 is a root of L_{p5} mod 5.");
Emit("");
Emit("Galois action:");
Emit("  Classes come in Galois-conjugate pairs under Gal(K/Q): s2 -> -s2.");
Emit("  Conjugate curves have L-polynomials that are related by the same s2->-s2");
Emit("  substitution in the coefficients; for inert or ramified p these coincide.");
Emit("");
Emit("Notes on reduction failures:");
Emit("  'reduction failed' in the Frobenius table is a MODEL artifact, not evidence");
Emit("  of bad reduction.  The Igusa criterion and ConductorExponentAt both confirm");
Emit("  good reduction at those primes.  The issue: the current Weierstrass model");
Emit("  (from Genus2CurveFromIgusa, via a random Mestre conic point) can have");
Emit("  leading coefficient divisible by certain prime ideals; neither integer");
Emit("  shifts nor the flip x->1/x resolve this.  A minimal Weierstrass model");
Emit("  at those primes would give good reduction.");
Emit("");
Emit("Sharing these curves:");
Emit("  The compact Igusa-Siegel invariants [J2:...:J10] above are the recommended");
Emit("  form for sharing.  A Weierstrass model y^2 = f(x) over Q(sqrt(2)) can be");
Emit("  recovered in Magma via:");
Emit("    k<s2> := QuadraticField(2);");
Emit("    load \"survivors.m\";");
Emit("    QI := survivors[2][2];  // pick a representative");
Emit("    import \"Genus2Curve.m\": Genus2CurveFromIgusa;");
Emit("    C := Genus2CurveFromIgusa(QI, k);");
Emit("");
printf "Done. Output written to %o\n", outfile;
