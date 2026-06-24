// factor_types.m -- Tabulate the factorization type of the Frobenius quartic h_P mod 11
// over all primes. Decides whether Jac(C)[11] is reducible 2+2 over K.
//   reducible 2+2  => h_P = (quadratic)(quadratic) ALWAYS: types in {[2,2],[1,1,2],[1,1,1,1]}
//                     (plus repeats). NEVER an irreducible quartic [4] or [1,3].
//   irreducible 4d (big image, ~Sp4) => [4] and [1,3] appear at a positive fraction.
//
// Run: magma factor_types.m

SetColumns(0);
load "mod11_frobdata.txt";
F11 := GF(11); R<x> := PolynomialRing(F11);

types := AssociativeArray();
n_4 := 0; n_13 := 0; examples4 := [];
for row in frobdata do
    cs := row[4];
    h := R ! [cs[5],cs[4],cs[3],cs[2],cs[1]];  // x^4 + c1 x^3 + c2 x^2 + c3 x + c4
    fac := Factorization(h);
    degs := Sort(&cat[ [Degree(t[1]) : j in [1..t[2]]] : t in fac ]);
    key := Sprintf("%o", degs);
    if not IsDefined(types, key) then types[key] := 0; end if;
    types[key] +:= 1;
    // flag irreducible quartic or [1,3] (impossible for a 2+2 rep)
    maxd := Max([Degree(t[1]) : t in fac]);
    if maxd ge 3 then
        if maxd eq 4 then n_4 +:= 1; else n_13 +:= 1; end if;
        if #examples4 lt 6 then Append(~examples4, <row[1],row[2],key>); end if;
    end if;
end for;

printf "Factorization types of h_P mod 11 over %o primes:\n", #frobdata;
for k in Sort([kk : kk in Keys(types)]) do
    printf "  %-14o : %o\n", k, types[k];
end for;
printf "\nirreducible-quartic [4] count : %o\n", n_4;
printf "type-[1,3] count              : %o\n", n_13;
if n_4 eq 0 and n_13 eq 0 then
    printf "=> NO degree>=3 irreducible factor ever. CONSISTENT with reducible 2+2.\n";
    printf "   (Jac(C)[11] splits as two 2-dim Galois reps over K.)\n";
else
    printf "=> degree>=3 factors appear (examples %o).\n", examples4;
    printf "   => Jac(C)[11] is NOT 2+2 reducible; the quartic does not give 2-dim subreps.\n";
end if;
