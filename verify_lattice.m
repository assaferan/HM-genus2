// verify_lattice.m -- check the hypothesis that lets Deligne-Serre lift our mod-l eigenvector
// to characteristic zero (results §4b, "From a mod-l eigenvector to a characteristic-0 form").
//
// The lifting lemma (Deligne-Serre 1974, Lemme 6.11) needs a finite-rank module over a DVR
// carrying commuting operators. Our space HilbertCuspForms(F, N, [2,2]) is already characteristic
// zero, so all that is required is a Hecke-STABLE lattice. Magma's Hecke matrices are rational and
// frequently NOT integral -- but if every denominator is prime to l, then the Z_(l)-span of Magma's
// basis is Hecke-stable and the lemma applies with R = Z_(l).
//
// That is exactly the same l-integrality that makes the reduction T_P -> T_P mod l well defined in
// kernel_torsion.m / grh_kernel.m, so this script also guards the arithmetic those rely on: were
// some denominator divisible by l, both the certificates AND the lifting argument would break.
//
// Usage:  magma verify_lattice.m              (default curves; seconds)
//         magma idx:=7 verify_lattice.m       (a single curve)
//         magma idx:=7 PBOUND:=60 verify_lattice.m
//
// Sentinel on success: "LATTICE: ALL PASS"
SetColumns(0);
load "torsion_data.m";

if not assigned PBOUND then PBOUND := 30; end if;
PBOUND := StringToInteger(Sprintf("%o", PBOUND));
// Defaults: one Q(sqrt2) curve and one Q(sqrt3) curve, both small. 4057.1 is deliberately
// included because ALL of its Hecke matrices are non-integral -- the interesting case.
if assigned idx then
    IDXS := [StringToInteger(Sprintf("%o", idx))];
else
    IDXS := [1, 28];
end if;

nbad := 0; nchecked := 0;
for i in IDXS do
    row := torsion_data[i];
    d := row[1]; l := row[2]; lab := row[3]; fc := row[4]; hc := row[5]; cc := row[7];
    K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K);
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

    M := HilbertCuspForms(K, base, [2,2]); dm := Dimension(M);
    nonint := 0; worst := 0; nP := 0;
    for pp in PrimesUpTo(PBOUND) do
        if pp eq 2 or pp eq l then continue; end if;
        for tup in Factorization(pp*OK) do
            P := tup[1];
            if Norm(base+P) ne 1 then continue; end if;
            T := Matrix(HeckeOperator(M, P)); nP +:= 1;
            D := LCM([ Denominator(z) : z in Eltseq(T) ]);
            if D ne 1 then nonint +:= 1; end if;
            v := Valuation(D, l);
            if v gt worst then worst := v; end if;
        end for;
    end for;

    nchecked +:= 1;
    if worst ne 0 then nbad +:= 1; end if;
    printf "%-10o d=%o l=%-3o dim=%-6o T_P checked=%-3o non-integral=%-3o max v_l(denom)=%o  %o\n",
        lab, d, l, dm, nP, nonint, worst,
        worst eq 0 select "l-integral" else "*** NOT l-integral ***";
end for;

printf "\n";
if nbad eq 0 then
    printf "Every Hecke matrix is l-integral, so the Z_(l)-span of Magma's basis is Hecke-stable\n";
    printf "and Deligne-Serre (1974, Lemme 6.11) applies with R = Z_(l).\n";
    printf "LATTICE: ALL PASS (%o curve(s))\n", nchecked;
else
    printf "LATTICE: FAIL (%o of %o curves have a denominator divisible by l)\n", nbad, nchecked;
end if;
exit;
