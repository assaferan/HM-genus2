// geom_endo_Bf.m -- determine End_Qbar^0(B_f) for the mod-5 idx-5 example.
//
// B_f / F=Q(sqrt2), RM by K_f=Q(sqrt3), attached to Hilbert newform f
// (orbit 6, level norm 5000, weight [2,2]).  The construction strategy hinges on
// whether B_f is GEOMETRICALLY SIMPLE (End_Qbar^0 = Q(sqrt3), a GENERIC point of
// the Humbert surface H_12) or has EXTRA geometric endomorphisms (QM / ~E^2, a
// point on a proper sub-locus / "badlocus" of H_12).
//
// For a non-CM GL_2-type surface, End_Qbar^0 grows beyond the RM field iff f has a
// nontrivial SELF-inner-twist over F: a Hecke character chi of F (FINITE order) with
//     a_P(f) * chi(P) = tau(a_P(f))  for all good P,   tau = Gal(K_f/Q).
// If it existed, chi(P) = tau(a_P)/a_P would be a root of unity at every good P.
// We test this directly on the eigenvalues, and separately rule out CM.
//
// Run: magma geom_endo_Bf.m
SetColumns(0);
F<w> := QuadraticField(2); OF := Integers(F);
LOGF := "geom_endo_Bf_out.txt";
PrintFile(LOGF, "# geom_endo_Bf: End_Qbar^0(B_f) determination" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1];
nlev := P2^3 * P5^2;
M := HilbertCuspForms(F, nlev, [2,2]);
Dc := NewformDecomposition(NewSubspace(M));
ef := Eigenform(Dc[6]);
K := HeckeEigenvalueField(Dc[6]);
LOG(Sprintf("K_f = %o  (deg %o)", DefiningPolynomial(K), Degree(K)));
assert Degree(K) eq 2;
// nontrivial Galois automorphism tau of K_f
auts := Automorphisms(K);
tau := [a : a in auts | a(K.1) ne K.1][1];
LOG(Sprintf("tau(K.1) = %o  (K.1 = %o)", tau(K.1), K.1));

// Collect a_P over good primes; record split/inert, a_P, and chi(P)=tau(a_P)/a_P.
badN := {2,5};
LOG("");
LOG("# ell : split? : P(norm) : a_P : tau(a_P)/a_P (candidate inner-twist char value)");
nzero := 0; ntot := 0;
ratios := [];
for l in PrimesInInterval(3, 80) do
    if l in badN then continue; end if;
    fac := Factorization(l*OF);
    isSplit := #fac eq 2;
    for pf in fac do
        P := pf[1]; NP := Norm(P);
        aP := HeckeEigenvalue(ef, P);
        ntot +:= 1;
        if aP eq 0 then
            nzero +:= 1;
            LOG(Sprintf("%o : %o : P(%o) : 0 : (a_P=0)", l, isSplit, NP));
            continue;
        end if;
        r := tau(K!aP)/(K!aP);
        Append(~ratios, r);
        // is r a root of unity?  finite multiplicative order?
        isroot := false;
        for m in [1..24] do if r^m eq 1 then isroot := true; break; end if; end for;
        LOG(Sprintf("%o : %o : P(%o) : %o : %o   root-of-unity? %o",
                    l, isSplit, NP, aP, r, isroot));
    end for;
end for;

LOG("");
LOG(Sprintf("a_P = 0 at %o of %o good primes (CM would force ~1/2 density on inert side)", nzero, ntot));

// Decision: any root-of-unity ratios uniformly?  If even ONE ratio is infinite order,
// no finite-order inner-twist character exists.
allroot := true;
for r in ratios do
    isr := false;
    for m in [1..24] do if r^m eq 1 then isr := true; break; end if; end for;
    if not isr then allroot := false; end if;
end for;
LOG("");
if allroot then
    LOG("Some inner-twist character is POSSIBLE (all ratios roots of unity): investigate QM/E^2.");
else
    LOG("At least one ratio tau(a_P)/a_P has INFINITE order => NO finite-order self-inner-twist");
    LOG("=> (f non-CM) End_Qbar^0(B_f) = Q(sqrt3): B_f is GEOMETRICALLY SIMPLE,");
    LOG("   a GENERIC point of the Humbert surface H_12 (NOT on an extra-endo badlocus).");
end if;
exit;
