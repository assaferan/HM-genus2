// debug_lambda.m -- diagnose why no consistent lambda was found.
SetColumns(0);
load "mod11_frobdata.txt";
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11); g11 := F11 ! PrimitiveRoot(11);
R11x<x> := PolynomialRing(F11);

function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;
function HbarFromCoeffs(cs) return R11x ! [cs[5],cs[4],cs[3],cs[2],cs[1]]; end function;
function Classify(hbar)
    fac := Factorization(hbar); lin := []; quad := R11x!0; hasq := false;
    for t in fac do
        if Degree(t[1]) eq 1 then
            r := -Coefficient(t[1],0)/Coefficient(t[1],1);
            for j in [1..t[2]] do Append(~lin, r); end for;
        elif Degree(t[1]) eq 2 then quad := t[1]; hasq := true; end if;
    end for;
    if hasq and #lin eq 2 then return "clean", lin, quad; end if;
    if #lin eq 4 then return "split", lin, R11x!0; end if;
    return "other", lin, R11x!0;
end function;

pbad := Factorization(66179*OK)[1][1];
p2   := Factorization(2*OK)[1][1];
M := (66179*OK)*(11*OK)*(3*OK)*(5*OK)*p2^8;
Rcl, mp := RayClassGroup(M, [1,2]);
invs := [ Order(Rcl.i) : i in [1..Ngens(Rcl)] ];   // ACTUAL generator orders (Eltseq basis)
printf "RCG invariants: %o,  |RCG| = %o\n", invs, &*invs;
printf "Ngens(Rcl) = %o,  [Order(Rcl.i)] = %o\n",
    Ngens(Rcl), [Order(Rcl.i) : i in [1..Ngens(Rcl)]];

// Construct chi_11 DIRECTLY: image of generator Rcl.i is N(mp(Rcl.i)) mod 11.
ng := Ngens(Rcl);
img_chi11 := [ F11 ! (Integers() ! (Norm(mp(Rcl.i)) mod 11)) : i in [1..ng] ];
printf "img_chi11 (from generator-ideal norms) = %o\n", img_chi11;
// test it against N(P) on clean primes
function evalg(img, es) return &*[ img[i]^es[i] : i in [1..ng] ]; end function;
bad := 0; tot := 0;
for row in frobdata do
    p := row[2]; rid := row[3];
    if p in {2,11} or p eq 66179 then continue; end if;
    typ := Classify(HbarFromCoeffs(row[4]));
    if typ ne "clean" then continue; end if;
    P := PrimeFromTag(p,rid); es := Eltseq(P @@ mp);
    tot +:= 1;
    if evalg(img_chi11, es) ne F11!(Integers()!(Norm(P) mod 11)) then bad +:= 1; end if;
end for;
printf "DIRECT chi_11 check: %o/%o clean primes MISMATCH (0 = machinery correct)\n\n", bad, tot;

// clean constraints
cons := [];
for row in frobdata do
    p := row[2]; rid := row[3];
    if p in {2,11} or p eq 66179 then continue; end if;
    typ, lin := Classify(HbarFromCoeffs(row[4]));
    if typ ne "clean" then continue; end if;
    P := PrimeFromTag(p, rid);
    Append(~cons, <Eltseq(P @@ mp), {lin[1],lin[2]}, F11!row[1], p>);
end for;
printf "clean constraints: %o\n", #cons;

function elsOrdDiv(t) return [ g11^((10 div t)*j) : j in [0..t-1] ]; end function;
choices := [ elsOrdDiv(Gcd(invs[i],10)) : i in [1..#invs] ];
CP := CartesianProduct(choices);
printf "number of order-|10 characters enumerated: %o\n\n", &*[#c : c in choices];

function evalchi(img, es)
    return &*[ img[i]^es[i] : i in [1..#invs] ];
end function;

// (a) SELF-TEST: can we reproduce chi_11 (value = N(P) mod 11)?  Search for it.
printf "SELF-TEST: searching for chi_11 (chi(P) = N(P) mod 11 on clean primes)...\n";
found_chi11 := false;
for tup in CP do
    img := [ tup[i] : i in [1..#invs] ];
    ok := true;
    for c in cons do
        if evalchi(img, c[1]) ne c[3] then ok := false; break; end if;  // c[3]=N(P) mod 11
    end for;
    if ok then found_chi11 := true; break; end if;
end for;
printf "  chi_11 reproducible by enumeration: %o\n\n", found_chi11;

// (b) best constraint-satisfaction for the lambda-pair condition
printf "Scanning all characters for best 'chi(P) in pair' satisfaction...\n";
best := -1; bestimg := [];
for tup in CP do
    img := [ tup[i] : i in [1..#invs] ];
    cnt := 0;
    for c in cons do
        if evalchi(img, c[1]) in c[2] then cnt +:= 1; end if;
    end for;
    if cnt gt best then best := cnt; bestimg := img; end if;
end for;
printf "  best satisfies %o / %o clean constraints\n", best, #cons;

// show where the best char fails (first few)
printf "\nFirst clean constraints vs best character:\n";
shown := 0;
for c in cons do
    v := evalchi(bestimg, c[1]);
    mark := v in c[2] select "ok" else "FAIL";
    printf "  p=%o N(P) mod11=%o pair=%o  chi=%o  [%o]\n", c[4], c[3], c[2], v, mark;
    shown +:= 1; if shown ge 12 then break; end if;
end for;
