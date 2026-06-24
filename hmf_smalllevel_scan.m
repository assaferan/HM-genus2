// hmf_smalllevel_scan.m -- Test the "level-lowered at 66179" hypothesis for the
// 2-dim mod-11 rep W of Curve A. W is unramified outside {2,11,66179}, weight (2,2),
// good at 11. If W is unramified at 66179 (level-lowering), its level is a pure power
// of the inert prime (2) (norm 4). Scan weight-(2,2) newforms over Q(sqrt5) of level
// (2)^e, e=1..6, and check each against the 132 exact a_p mod 11 in mod11_apdata.txt.
//
// Match rule: for each rational prime p (good, <= PBound) in our data, the multiset
// { a_P mod 11 : P | p } from W must be contained in the form's { a_P mod 11 : P | p }.
// Tolerant of a global sign flip (Frobenius-trace convention).
//
// Run: magma hmf_smalllevel_scan.m

K<w> := QuadraticField(5);
OK := Integers(K);
F11 := GF(11);
PBound := 200;

// ---- load W data: aplist := [ <N(p), p, a_p mod 11>, ... ] ----
load "mod11_apdata.txt";
// group W values by rational prime p
Wbyp := AssociativeArray();
for t in aplist do
    p := t[2];
    if not IsDefined(Wbyp, p) then Wbyp[p] := []; end if;
    Append(~Wbyp[p], F11 ! t[3]);
end for;
Wprimes := Sort([ p : p in Keys(Wbyp) | p le PBound and p ne 11 ]);
printf "Loaded W data: %o rational primes p <= %o with exact a_p mod 11.\n\n",
    #Wprimes, PBound;

p2 := Factorization(2*OK)[1][1];   // the inert prime (2), norm 4

// Reduce a Hecke eigenvalue (rational or in a number field Hecke field) to F_11.
// Returns <true, value> if there is a degree-1 prime above 11 in the Hecke field
// (required for an F_11-valued match), else <false, _>.
function RedMod11(a)
    if a in Rationals() then return true, F11 ! a; end if;
    Kf := Parent(a); OKf := MaximalOrder(Kf);
    for fac in Factorization(11*OKf) do
        P := fac[1];
        Ff, red := ResidueClassField(P);
        if #Ff eq 11 then
            // identify Ff with F_11
            return true, F11 ! (Integers() ! red(OKf ! a));
        end if;
    end for;
    return false, _;
end function;

function FormApMultiset(eig, p)
    // multiset of a_P mod 11 over primes P | p; second return false if not F_11-reducible
    res := [];
    for fac in Factorization(p*OK) do
        ok, v := RedMod11(HeckeEigenvalue(eig, fac[1]));
        if not ok then return res, false; end if;
        Append(~res, v);
    end for;
    return res, true;
end function;

function ContainedAsMultiset(small, big)   // is multiset small subset of big?
    b := big;
    for s in small do
        idx := Index(b, s);
        if idx eq 0 then return false; end if;
        Remove(~b, idx);
    end for;
    return true;
end function;

for e in [1..6] do
    N := p2^e;
    printf "=== level (2)^%o, norm %o ===\n", e, Norm(N);
    M := HilbertCuspForms(K, N, [2,2]);
    Mnew := NewSubspace(M);
    dnew := Dimension(Mnew);
    printf "  new subspace dimension = %o\n", dnew;
    if dnew eq 0 then printf "  (no newforms)\n\n"; continue; end if;
    decomp := NewformDecomposition(Mnew);
    printf "  %o Galois orbits of newforms\n", #decomp;
    for nf in decomp do
        eig := Eigenform(nf);
        // only rational eigenform orbits give a chance at an F_11 match here, but
        // check all: reduce eigenvalues mod a prime above 11 in the Hecke field.
        ok := true; ok_neg := true; checked := 0; reducible := true;
        for p in Wprimes do
            if p eq 2 then continue; end if;
            formset, red := FormApMultiset(eig, p);
            if not red then reducible := false; break; end if;
            Wset := Wbyp[p];
            if not ContainedAsMultiset(Wset, formset) then ok := false; end if;
            if not ContainedAsMultiset([-a : a in Wset], formset) then ok_neg := false; end if;
            checked +:= 1;
            if not ok and not ok_neg then break; end if;
        end for;
        if not reducible then
            printf "    orbit dim %o: no degree-1 prime above 11 in Hecke field -> cannot match F_11 rep\n",
                Dimension(nf);
            continue;
        end if;
        hf := HeckeEigenvalueField(nf);
        printf "    orbit dim %o, Hecke field deg %o: match=%o  match(-)=%o  (checked %o primes)\n",
            Dimension(nf), Degree(hf), ok, ok_neg, checked;
        if ok or ok_neg then
            printf "    *** CANDIDATE MATCH at level norm %o ***\n", Norm(N);
        end if;
    end for;
    printf "\n";
end for;
