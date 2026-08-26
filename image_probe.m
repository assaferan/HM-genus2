// image_probe.m -- cheap statistical probe of the image of sigma, to see whether any example
// admits an UNCONDITIONAL modularity route (and hence a certificate free of GRH).
//
// Why the trace-zero frequency is the discriminator. sigma is 2-dimensional with det = chi_l.
//   * If sigma is DIHEDRAL (induced from a character of a quadratic extension K/F), then
//     tr sigma(Frob_P) = 0 for EXACTLY the primes P inert in K -- density 1/2. Such a sigma is
//     automorphic by automorphic induction (a theta series), with NO conjecture: modularity
//     becomes unconditional, and identifying which form then needs only a Sturm bound.
//   * If the image contains SL_2(F_l), traces are roughly equidistributed over F_l, so
//     tr = 0 has density about 1/l -- 9% at l=11, 8% at l=13.
// The two predictions differ by a factor of ~6, so a few hundred primes settle it decisively.
//
// We use the fingerprint we already trust: tr sigma(Frob_P) = -#C(F_P) (mod l)  (§1).
//
// Usage:  magma image_probe.m              (all 39 curves)
//         magma idx:=3 BOUND:=2000 image_probe.m
//
// Sentinel: "IMAGE PROBE: DONE"
SetColumns(0);
load "torsion_data.m";

if not assigned BOUND then BOUND := 1000; end if;
BOUND := StringToInteger(Sprintf("%o", BOUND));
if assigned idx then
    IDXS := [StringToInteger(Sprintf("%o", idx))];
else
    IDXS := [1..#torsion_data];
end if;

printf "trace-zero frequency of sigma  (good primes N(P) <= %o)\n", BOUND;
printf "dihedral predicts ~0.500; image containing SL_2(F_l) predicts ~1/l\n\n";
printf "%-11o %-3o %-4o %6o %7o %8o %8o  %o\n",
       "label", "l", "1/l", "primes", "zeros", "observed", "distinct", "reading";

ndih := 0;
for i in IDXS do
    row := torsion_data[i];
    d := row[1]; l := row[2]; lab := row[3]; fc := row[4]; hc := row[5]; cc := row[7];
    K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K); Fl := GF(l);
    f := R![K| c[1]+c[2]*a : c in fc];
    hh := R![K| c[1]+c[2]*a : c in hc];
    C := HyperellipticCurve(f, hh);
    if d eq 2 then
        base := ideal<OK | cc[1]+cc[2]*a>;
    else
        base := 1*OK;
        for pe in Factorization(Conductor(C)) do
            if Norm(pe[1]) mod 2 ne 0 then base := base*pe[1]^pe[2]; end if;
        end for;
    end if;

    nP := 0; nz := 0; vals := {};
    for pp in PrimesUpTo(BOUND) do
        if pp eq 2 or pp eq l then continue; end if;
        for tup in Factorization(pp*OK) do
            P := tup[1];
            if Norm(P) gt BOUND or Norm(base+P) ne 1 then continue; end if;
            kk, red := ResidueClassField(P); Rk := PolynomialRing(kk);
            ok := true; Ck := 0;
            try
                Ck := HyperellipticCurve(Rk![red(c) : c in Coefficients(f)],
                                          Rk![red(c) : c in Coefficients(hh)]);
            catch e
                ok := false;
            end try;
            if not ok then continue; end if;
            t := Fl!(-#Ck);
            nP +:= 1; Include(~vals, t);
            if t eq 0 then nz +:= 1; end if;
        end for;
    end for;

    frac := nP eq 0 select 0.0 else RealField(4)!(nz/nP);
    expected := RealField(4)!(1/l);
    // Call it dihedral only if the zero-frequency is far above 1/l and near 1/2.
    isdih := (nP ge 50) and (frac gt 0.30);
    if isdih then ndih +:= 1; end if;
    printf "%-11o %-3o %-4.3o %6o %7o %8.3o %8o  %o\n",
        lab, l, expected, nP, nz, frac, #vals,
        isdih select "*** DIHEDRAL CANDIDATE ***"
               else (#vals eq l select "all l traces occur -- large image" else "consistent with large image");
end for;

printf "\n";
if ndih eq 0 then
    printf "No dihedral candidate: every sigma has trace-zero frequency near 1/l, not 1/2.\n";
    printf "So none of these examples is automorphic by induction, and the easy unconditional\n";
    printf "route is unavailable -- the GRH-conditional certificate stands as the statement.\n";
else
    printf "%o DIHEDRAL CANDIDATE(S) -- these may be unconditionally modular (theta series);\n", ndih;
    printf "confirm by exhibiting the quadratic character and then use a Sturm bound.\n";
end if;
printf "IMAGE PROBE: DONE\n";
exit;
