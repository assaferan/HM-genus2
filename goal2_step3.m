// goal2_step3.m -- Goal 2, Step 3: pin down the images of sigma_A and rho_f in
// GL_2(F_49) precisely (needed to structure the Faltings-Serre certificate).
//
// Dickson: an irreducible subgroup of GL_2(F_{7^f}) containing a transvection
// (element of order 7) contains SL_2(F_{7^m}), where F_{7^m} = field generated
// by the TRACES of the group.  We show the trace field is F_49 for BOTH reps,
// hence both images contain SL_2(F_49).  We also census Frobenius element types
// (central / unipotent / split-ss / nonsplit-ss) -- the data for covering the
// conjugacy classes in the FS test set.
//
// Run: magma goal2_step3.m

SetColumns(0);
F<w> := QuadraticField(10); OF := Integers(F);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
fpoly := -10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9;
C := HyperellipticCurve(fpoly);
F7 := GF(7); F49<z> := GF(49); PZ<Z> := PolynomialRing(F49);
badp := {2,3,5,7};

LOGF := "goal2_step3_out.txt";
PrintFile(LOGF, "# goal2_step3: precise images of sigma_A, rho_f in GL_2(F_49)" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

// classify (tr,det) in F49: element type of a matrix with this char poly
function EltType(tr, det)
    disc := tr^2 - 4*det;
    if disc eq 0 then
        return (tr eq 0) select "central" else "unipotent(x scalar)";
    elif IsSquare(disc) then return "split-ss";
    else return "nonsplit-ss"; end if;
end function;

// ---------- sigma_A: traces from A's L-polynomials ----------
LOG("=== sigma_A (from A[7]) ===");
sigA := [];    // <tr in F49, det in F49>
for l in PrimesInInterval(3, 600) do
    if l in badp then continue; end if;
    LA := LPolynomial(ChangeRing(C, GF(l)));
    a1 := -Coefficient(LA,1); a2 := Coefficient(LA,2);
    if IsSquare(GF(l)!10) then          // split: two sigma-traces, det = l
        for rt in Roots(Z^2 - F49!(F7!a1)*Z + (F49!(F7!a2) - 2*F49!(F7!l))) do
            for k in [1..rt[2]] do Append(~sigA, <rt[1], F49!(F7!l)>); end for;
        end for;
    else                                 // inert: tr = -a2, det = l^2
        Append(~sigA, < F49!(F7!(-a2)), F49!(F7!(l^2)) >);
    end if;
end for;
trfieldA := sub< F49 | [e[1] : e in sigA] >;
has7A := exists{ e : e in sigA | e[1]^2 eq 4*e[2] and e[1] ne 0 };
LOG(Sprintf("  %o Frobenius traces; trace field = F_%o", #sigA, #trfieldA));
LOG(Sprintf("  transvection (order-7) present: %o", has7A));
cnt := AssociativeArray();
for e in sigA do t := EltType(e[1],e[2]); cnt[t] := (IsDefined(cnt,t) select cnt[t] else 0)+1; end for;
LOG(Sprintf("  element-type census: %o", [<k, cnt[k]> : k in Keys(cnt)]));

// ---------- rho_f: traces from the eigenform (split primes: single-valued) ----------
LOG("\n=== rho_f (orbit 10 mod lambda|7) ===");
P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1]; P3a := Factorization(3*OF)[1][1];
M := HilbertCuspForms(F, P2^2*P3a*P5^3, [2,2]);
D := NewformDecomposition(NewSubspace(M));
e := Eigenform(D[10]); E := HeckeEigenvalueField(D[10]); OE := Integers(E);
found := false;
for pl in Factorization(7*OE) do
    Fq0, r0 := ResidueClassField(pl[1]);
    if Degree(Fq0) eq 2 then red := r0; found := true; break; end if;
end for;
assert found;
// (only the trace FIELD of rho_f is needed here -- whether traces leave F_7 is
//  basis-independent, so no cross-field embedding is required.)
rhof := [];    // <tr in F49, det=Norm(P) in F49>
for l in PrimesInInterval(3, 600) do
    if l in badp or not IsSquare(GF(l)!10) then continue; end if;   // split primes: norm l
    for P in [pl[1] : pl in Factorization(l*OF)] do
        a := red(OE!HeckeEigenvalue(e, P));   // in the eigenform's F_49
        // move to our F49 by matching: a is in GF(49); coerce via its min poly over F7
        tr := F49 ! Eltseq(a, F7);            // same F7-basis coordinates
        Append(~rhof, < tr, F49!(F7!l) >);
    end for;
end for;
trfieldF := sub< F49 | [e2[1] : e2 in rhof] >;
LOG(Sprintf("  %o traces (split primes); trace field = F_%o", #rhof, #trfieldF));

LOG("");
LOG(Sprintf("CONCLUSION: sigma_A trace field F_%o, rho_f trace field F_%o.", #trfieldA, #trfieldF));
LOG("If both = F_49 and both irreducible with a transvection => both images");
LOG("contain SL_2(F_49) (Dickson).  det = chi_7 surjects onto F_7^*, so");
LOG("im = { g in GL_2(F_49) : det g in F_7^* } (mod exact index) -- BIG, matching.");
exit;
