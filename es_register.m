// es_register.m -- compute the data needed to register B_f (2.2.8.1-5000.1-j) in
// EichlerShimuraHMF/src/Labels.m: level generators, dimension, Hecke field, and a trace
// fingerprint (LMFDB ideal label -> Trace(a_P)) that UNIQUELY picks our orbit among the
// dim-2 orbits.  Uses the package's own field/ideal conventions for consistency.
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
AttachSpec("../EichlerShimuraHMF/src/spec");

LOGF := "es_register_out.txt";
PrintFile(LOGF, "# registration data for 2.2.8.1-5000.1-j" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

F<w> := LMFDBField("2.2.8.1");   // = Q(sqrt2), w = sqrt2
OF := Integers(F);
LOG(Sprintf("F = %o, w^2 = %o", DefiningPolynomial(F), w^2));

// level = P2^3 * P5^2 (norm 5000); 2 ramifies (P2=(w)), 5 is inert (P5=(5))
P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1];
NN := P2^3 * P5^2;
LOG(Sprintf("level norm %o, principal? %o", Norm(NN), IsPrincipal(NN)));
b, gen := IsPrincipal(NN);
if b then LOG(Sprintf("  generator = %o = %o", gen, Eltseq(gen))); end if;
// cross-check ideal from [[0,50]] (=50*w)
NN2 := ideal<OF | [F!50*w]>;
LOG(Sprintf("ideal(50w) == level? %o  (norm %o)", NN2 eq NN, Norm(NN2)));

// build space, decompose, find dim-2 orbits
M := HilbertCuspForms(F, NN);
D := NewformDecomposition(NewSubspace(M));
LOG(Sprintf("#orbits = %o", #D));
dim2 := [elt : elt in D | Dimension(elt) eq 2];
LOG(Sprintf("#dim-2 orbits = %o", #dim2));
fs := [* Eigenform(elt) : elt in dim2 *];
for i->elt in dim2 do
    Ki := HeckeEigenvalueField(elt);
    LOG(Sprintf("  dim-2 orbit %o: Hecke field %o", i, DefiningPolynomial(Ki)));
end for;
// our B_f: Hecke field x^2 + 2x - 2
target := PolynomialRing(Rationals()).1^2 + 2*PolynomialRing(Rationals()).1 - 2;
ours := [i : i->elt in dim2 | DefiningPolynomial(HeckeEigenvalueField(elt)) eq target];
LOG(Sprintf("orbit(s) with Hecke field x^2+2x-2: %o", ours));

// trace fingerprint at small primes: LMFDB label -> Trace(a_P), for EACH dim-2 orbit,
// to find a discriminating set that uniquely picks ours.
LOG("");
LOG("# trace table: prime label (norm) : [Trace(a_P) for each dim-2 orbit]");
smallprimes := [];
for l in [3,7,9,17,23,25,31] do
    for pf in Factorization(l*OF) do P := pf[1];
        if Norm(P) ne l then continue; end if;
        Append(~smallprimes, P);
    end for;
end for;
for P in smallprimes do
    lab := LMFDBLabel(P);
    traces := [ Trace(HeckeEigenvalue(f, P)) : f in fs ];
    LOG(Sprintf("  %o (norm %o) : %o", lab, Norm(P), traces));
end for;
exit;
