// goal2_innertwist.m -- inner-twist structure of the degree-18 Hilbert newform f
// (orbit 10, F=Q(sqrt10), level norm 1500, weight (2,2)).  For a totally real f
// with trivial nebentypus an inner twist is (gamma, chi) with gamma in Aut(E) and
// chi a QUADRATIC Hecke character such that gamma(a_P) = chi(P) a_P for all P.
// The fixed field E_0 of the inner-twist group Gamma gives the effective building
// block: B_f ~ B_0^{|Gamma|} with dim B_0 = [E_0:Q].  If [E_0:Q] is small, the
// unconditional torsion-field comparison (F_f side) becomes tractable.
//
// Run: magma goal2_innertwist.m
SetColumns(0);
F<w> := QuadraticField(10); OF := Integers(F);
LOGF := "goal2_innertwist_out.txt";
PrintFile(LOGF, "# inner twists of orbit-10 (deg-18 Hilbert newform)" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1]; P3a := Factorization(3*OF)[1][1];
level := P2^2 * P3a * P5^3;
LOG("building HilbertCuspForms(1500,[2,2]) ...");
t0 := Cputime();
M := HilbertCuspForms(F, level, [2,2]);
D := NewformDecomposition(NewSubspace(M));
e := Eigenform(D[10]); E := HeckeEigenvalueField(D[10]);
assert Degree(E) eq 18;
LOG(Sprintf("  orbit 10, Hecke field E of degree %o (%.1o s)", Degree(E), Cputime(t0)));
LOG(Sprintf("  E totally real: %o", IsTotallyReal(E)));

// automorphisms of E
auts := Automorphisms(E);
LOG(Sprintf("  #Aut(E) = %o", #auts));

// Hecke eigenvalues a_P in E at split good primes (norm l, cheap)
aps := [];   // <P, a_P>
for l in PrimesInInterval(3, 250) do
    if l in {2,5,7} then continue; end if;
    if not IsSquare(GF(l)!10) then continue; end if;    // split -> norm l
    for P in [pl[1] : pl in Factorization(l*OF)] do
        Append(~aps, < P, HeckeEigenvalue(e, P) >);
    end for;
end for;
LOG(Sprintf("  computed %o Hecke eigenvalues a_P (split primes)", #aps));

// for each automorphism gamma, test whether gamma(a_P) = eps_P * a_P with eps_P in {+-1}
LOG("\n  automorphism | inner twist? | sign pattern (on nonzero a_P)");
innergammas := [* *];
for gi -> g in auts do
    signs := []; ok := true; nz := 0;
    for pr in aps do
        a := pr[2];
        if a eq 0 then continue; end if;
        nz +:= 1;
        r := g(a)/a;                 // should be +-1 for an inner twist
        if r eq 1 then Append(~signs, 1);
        elif r eq -1 then Append(~signs, -1);
        else ok := false; break; end if;
    end for;
    if ok and nz ge 8 then
        Append(~innergammas, g);
        np := #[s : s in signs | s eq 1]; nm := #[s : s in signs | s eq -1];
        LOG(Sprintf("  aut %2o        | YES          | +1 at %o, -1 at %o (of %o)", gi, np, nm, nz));
    end if;
end for;
LOG(Sprintf("\n  inner-twist group order |Gamma| = %o", #innergammas));

// fixed field E_0 of the inner-twist automorphisms
if #innergammas gt 0 then
    E0 := FixedField(E, [g : g in innergammas]);
    LOG(Sprintf("  fixed field E_0 = E^Gamma has degree [E_0:Q] = %o", Degree(E0)));
    LOG(Sprintf("  => B_f ~ B_0^{%o} with dim B_0 = %o", #innergammas, Degree(E0)));
    if Degree(E0) le 3 then
        LOG("  ==> SMALL building block: the F_f torsion-field computation is tractable.");
    else
        LOG("  ==> building block still of dimension " cat IntegerToString(Degree(E0)) cat ".");
    end if;
    LOG(Sprintf("  E_0 defining polynomial: %o", DefiningPolynomial(E0)));
end if;
LOG("DONE.");
exit;
