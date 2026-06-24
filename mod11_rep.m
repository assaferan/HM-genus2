// mod11_rep.m -- Extract the 2-dimensional mod-11 Galois representation W = V^perp/V
// from Curve A (the Q(sqrt(5)) (1,11) Gross-Popescu curve), and tabulate its
// Frobenius traces a_p mod 11.
//
// Background. A = Jac(Curve A) has a rational cyclic 11-isogeny, i.e. a Galois-stable
// isotropic line V subset A[11] (a character lambda). Under the Weil pairing
// V subset V^perp (dim 3), and W := V^perp/V is a 2-dim F_11-rep with det W = chi_11.
// Good reduction at both primes above 11 => W is Barsotti-Tate there => Serre weight (2,2).
//
// Extraction at a good prime p (p != 11, p not bad):
//   h_p(x) = char poly of Frob_p on the 4-dim Tate module (deg 4, reverse of L_p).
//   Mod 11, h_p = (x - lambda)(x - N(p)/lambda) * (2-dim factor).
//   Since V is rational, lambda, lambda' = N(p)/lambda are in F_11 (>=2 linear factors).
//   If W is irreducible at p: the remaining factor is an irreducible quadratic
//     x^2 - a_p x + N(p)  (constant term == N(p) mod 11)  ->  read a_p mod 11.
//   If W splits at p: h_p is totally split; flagged AMBIGUOUS (needs global lambda).
//
// Output: table (N(p), p, a_p mod 11, factor type) and a_p list to mod11_apdata.txt.
//
// Run: magma mod11_rep.m

SetColumns(0);   // disable line-wrapping so PrintFile output stays Magma-loadable
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, ConductorExponentAt, BadPrimesFromInvariants, PrimesAbove;

NormBound := 2000;   // include good primes p with N(p) <= NormBound

k<s5> := QuadraticField(5);
OK := Integers(k);
ZT<T> := PolynomialRing(Integers());
F11 := GF(11);
R11<x> := PolynomialRing(F11);

// ---- Curve A Igusa invariants (= survivors 1,2 of survivors_qsqrt5_gp11.m) ----
QI := [k|
    1,
    1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895)
];

printf "Reconstructing Curve A over K = Q(sqrt(5))...\n";
C := Genus2CurveFromIgusa(QI, k);
Cmin, d := MinimalTwist(C, k);
printf "  minimal twist factor d = %o\n", d;
bad := BadPrimesFromInvariants(QI, k);
bad_set := {p : p in bad};
printf "  bad primes (odd part): %o\n", [Norm(p) : p in bad];
printf "  11 splits: primes above 11 have norms %o\n\n",
    [Norm(p) : p in PrimesAbove(k, 11)];

// ---- integral model + mod-P reduction (same as verify_curveA.m) ----
function MakeIntegral(CC)
    f0, h0 := HyperellipticPolynomials(CC);
    K0 := BaseRing(CC);
    m := LCM([Integers()| Denominator(c) : c in Coefficients(f0) cat Coefficients(h0) cat [K0!1]]);
    return HyperellipticCurve(Parent(f0)!(m^2*f0), Parent(f0)!(m*h0));
end function;

function TryReduceModP(CI, P)
    K0 := BaseRing(CI); OK0 := Integers(K0);
    Fp, redP := ResidueClassField(P); Rx<xx> := PolynomialRing(Fp);
    f0, _ := HyperellipticPolynomials(CI); n := Degree(f0); R := Parent(f0);
    for flip in [0, 1] do
        f := flip eq 0 select f0
             else Polynomial([Coefficient(f0, n-i) : i in [0..n]]);
        for sh in [0..20] do
            fp := Rx![redP(OK0!c) : c in Coefficients(Evaluate(f, R.1 + sh))];
            if Degree(fp) eq n then return true, HyperellipticCurve(fp, Rx!0); end if;
        end for;
    end for;
    return false, _;
end function;

CI := MakeIntegral(Cmin);

// ---- main loop: extract a_p mod 11 of the 2-dim factor W ----
// Frobenius char poly h_p(x) = sum c_i x^(4-i) where L_p = sum c_i T^i (reverse).
function FrobCharPoly(lp)   // lp = L_p(T) in Z[T], degree 4
    cs := Coefficients(lp);            // [c0..c4], c0 = 1
    return ZT ! [cs[5], cs[4], cs[3], cs[2], cs[1]];  // x^4 + c1 x^3 + ... + c4
end function;

printf "%-7o %-4o %-6o %-7o %o\n", "N(p)", "p", "Np%11", "a_p%11", "factor type mod 11";
printf "%o\n", "-"^70;

aplist := [];   // <Norm(p), p, a_p mod 11>  (only unambiguous)
ambiguous := [];
frobdata := [];  // <Norm(p), p, [c0..c4]>  raw L_p(T) coeffs -- reusable for any mod-ell work
pp := 2;
while pp le NormBound do
    if pp eq 11 then pp := NextPrime(pp); continue; end if;
    for P in PrimesAbove(k, pp) do
        if Norm(P) gt NormBound then continue; end if;
        if P in bad_set then continue; end if;
        if ConductorExponentAt(Cmin, P) ne 0 then continue; end if;
        ok, Cp := TryReduceModP(CI, P);
        if not ok then continue; end if;
        Np := Norm(P);
        // residue identifier: for a degree-1 prime (Np = pp), store red(s5) in 0..pp-1
        // to distinguish the two primes above a split p; for inert (Np = pp^2), store -1.
        if Np eq pp then
            Fp, redp := ResidueClassField(P);
            rid := Integers() ! (Fp ! redp(OK ! s5));
        else
            rid := -1;
        end if;
        lp := ZT ! LPolynomial(Cp);
        Append(~frobdata, <Np, pp, rid, Coefficients(lp)>);
        h  := FrobCharPoly(lp);
        hbar := R11 ! h;
        fac := Factorization(hbar);
        // describe factor-degree multiset
        degs := &cat[ [Degree(t[1]) : j in [1..t[2]]] : t in fac ];
        Sort(~degs);
        // find an irreducible quadratic factor with constant term == N(p) mod 11
        Npbar := F11 ! Np;
        found := false; ap := F11!0; q := hbar;
        for t in fac do
            qq := t[1];
            if Degree(qq) eq 2 and IsIrreducible(qq) and Coefficient(qq,0) eq Npbar then
                // qq = x^2 + b x + c  ==  x^2 - a_p x + N(p)  => a_p = -b
                ap := -Coefficient(qq, 1);
                found := true; q := qq;
                break;
            end if;
        end for;
        ftype := Sprintf("%o", degs);
        if found then
            printf "%-7o %-4o %-6o %-7o %o  W=(%o)\n",
                Np, pp, Integers()!Npbar, Integers()!ap, ftype, q;
            Append(~aplist, <Np, pp, Integers()!ap>);
        else
            printf "%-7o %-4o %-6o %-7o %o  [AMBIGUOUS: W split or repeated]\n",
                Np, pp, Integers()!Npbar, "?", ftype;
            Append(~ambiguous, <Np, pp, hbar>);
        end if;
    end for;
    pp := NextPrime(pp);
end while;

printf "\n%o unambiguous a_p, %o ambiguous (totally split) primes.\n",
    #aplist, #ambiguous;

// ---- write a_p data for the HMF matching step ----
fname := "mod11_apdata.txt";
PrintFile(fname, "// Curve A: 2-dim mod-11 rep W = V^perp/V over Q(sqrt(5))" : Overwrite := true);
PrintFile(fname, "// rows: <N(p), rational p, a_p mod 11>  (W irreducible at p)");
PrintFile(fname, Sprintf("aplist := %o;", aplist));
printf "\nWrote %o (%o entries).\n", fname, #aplist;

// ---- write raw L-polynomial data (reusable) ----
gname := "mod11_frobdata.txt";
PrintFile(gname, "// Curve A: raw Frobenius L-polynomials over Q(sqrt(5)), good primes p != 11" : Overwrite := true);
PrintFile(gname, "// rows: <N(p), rational p, rid, [c0,c1,c2,c3,c4]>;  L_p(T)=sum c_i T^i.");
PrintFile(gname, "// rid = red(sqrt5) mod p in 0..p-1 identifies the prime above split p; -1 if inert.");
PrintFile(gname, Sprintf("frobdata := %o;", frobdata));
printf "Wrote %o (%o entries).\n", gname, #frobdata;
