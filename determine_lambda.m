// determine_lambda.m -- Determine the isogeny character lambda of Curve A's
// 2-dim mod-11 rep W, then use it to resolve the totally-split primes.
//
// lambda: Gal(K-bar/K) -> F_11^* is the Galois action on the rational line V (the
// 11-isogeny kernel). It is a finite-order (order | 10) Hecke character of K=Q(sqrt5),
// unramified outside {2,11,66179}. At a "clean" prime P (W irreducible), h_P mod 11 has
// exactly two linear roots = {lambda(Frob_P), lambda'(Frob_P)} (lambda' = N(P)/lambda).
// We find lambda as a hom RayClassGroup -> F_11^* consistent with all clean pairs, then
// at each totally-split prime remove the lambda-pair from the 4 roots to read off W.
//
// Run: magma determine_lambda.m   (needs mod11_frobdata.txt from mod11_rep.m)

SetColumns(0);   // keep PrintFile output Magma-loadable
load "mod11_frobdata.txt";   // frobdata := [ <Np, p, rid, [c0..c4]>, ... ]
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
g11 := F11 ! PrimitiveRoot(11);   // generator of F_11^* (=2)
R11x<x> := PolynomialRing(F11);

// reconstruct the prime ideal P from <p, rid>
function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;     // inert
    return ideal<OK | p, w - rid>;                       // degree-1 prime
end function;

// h_P(x) mod 11 from L_p coeffs [c0..c4]: h = x^4 + c1 x^3 + c2 x^2 + c3 x + c4
function HbarFromCoeffs(cs)
    return R11x ! [cs[5], cs[4], cs[3], cs[2], cs[1]];
end function;

// classify a prime's mod-11 Frobenius poly. Returns <type, linroots, Wquad-or-0>
// type: "clean" (irred quadratic factor present), "split" (>=3 linear), "other".
function Classify(hbar)
    fac := Factorization(hbar);
    linroots := [];  quad := R11x!0; hasquad := false;
    for t in fac do
        if Degree(t[1]) eq 1 then
            r := -Coefficient(t[1],0)/Coefficient(t[1],1);
            for j in [1..t[2]] do Append(~linroots, r); end for;
        elif Degree(t[1]) eq 2 and IsIrreducible(t[1]) then
            quad := t[1]; hasquad := true;
        end if;
    end for;
    if hasquad and #linroots eq 2 then return "clean", linroots, quad; end if;
    if #linroots eq 4 then return "split", linroots, R11x!0; end if;
    return "other", linroots, R11x!0;
end function;

// ---- gather clean-prime constraints: <class in R, {pair}> ----
pbad := Factorization(66179*OK)[1][1];
p2   := Factorization(2*OK)[1][1];

function FindLambda(mfin)
    Rcl, mp := RayClassGroup(mfin, [1,2]);
    // NB: use the orders of the ACTUAL generators Rcl.i (these are what Eltseq uses),
    // NOT Invariants(Rcl) -- they generally differ.
    ng := Ngens(Rcl);
    ords := [ Order(Rcl.i) : i in [1..ng] ];
    // build clean constraints
    cons := [];
    for row in frobdata do
        p := row[2]; rid := row[3];
        if p in {2,11} or (p eq 66179) then continue; end if;
        typ, lin, _ := Classify(HbarFromCoeffs(row[4]));
        if typ ne "clean" then continue; end if;
        P := PrimeFromTag(p, rid);
        cl := P @@ mp;
        Append(~cons, <Eltseq(cl), {lin[1], lin[2]}>);
    end for;
    // enumerate homs Rcl -> F_11^* of order dividing 10
    function elsOrdDiv(t) return [ g11^((10 div t)*j) : j in [0..t-1] ]; end function;
    choices := [ elsOrdDiv(Gcd(ords[i],10)) : i in [1..ng] ];
    CP := CartesianProduct(choices);
    good := [];
    for tup in CP do
        img := [ tup[i] : i in [1..ng] ];
        okc := true;
        for c in cons do
            es := c[1];
            val := ng eq 0 select F11!1 else &*[ img[i]^es[i] : i in [1..ng] ];
            if val notin c[2] then okc := false; break; end if;
        end for;
        if okc then Append(~good, img); end if;
    end for;
    return good, Rcl, mp, ords, #cons;
end function;

// use FULL ideals; include 3,5 too (Mestre model "fails" at 3,5 -- A may be genuinely
// bad there, so lambda could ramify at those primes). 5 ramifies in K.
base := (66179*OK)*(11*OK)*(3*OK)*(5*OK);
good := []; mfin := base;
for k in [4,6,8,10,12] do
    mfin := base * p2^k;
    printf "Searching for lambda over modulus pbad * (11) * p2^%o ...\n", k;
    good, Rcl, mp, invs, ncons := FindLambda(mfin);
    printf "  RCG invariants %o, %o constraints, consistent chars: %o\n", invs, ncons, #good;
    if #good gt 0 then break; end if;
end for;

if #good eq 0 then
    printf "  NO consistent lambda even at p2^12 -- reconsider assumptions. Stop.\n";
    exit;
end if;

// pick lambda = first consistent character; define evaluator
lam_img := good[1];
function chi(img, P)
    cl := Eltseq(P @@ mp);
    return #invs eq 0 select F11!1 else &*[ img[i]^cl[i] : i in [1..#invs] ];
end function;

// order of lambda
ord := 1;
for a in lam_img do ord := Lcm(ord, Order(a)); end for;
printf "  lambda has order %o.  (lambda' = chi_11/lambda is also consistent: pair-symmetry)\n", ord;

// sanity: lambda(Frob_P) lands in the clean pair, AND lambda*lambda' = N(P)
printf "\nSanity check on clean primes (lambda in pair, lambda*lambda' = N(P) mod 11):\n";
nbad := 0; ncheck := 0;
for row in frobdata do
    p := row[2]; rid := row[3];
    if p in {2,11} or p eq 66179 then continue; end if;
    typ, lin, _ := Classify(HbarFromCoeffs(row[4]));
    if typ ne "clean" then continue; end if;
    P := PrimeFromTag(p, rid);
    lv := chi(lam_img, P);
    Npbar := F11 ! row[1];
    ncheck +:= 1;
    if lv notin {lin[1],lin[2]} or lin[1]*lin[2] ne Npbar then nbad +:= 1; end if;
end for;
printf "  checked %o clean primes, inconsistencies: %o\n", ncheck, nbad;

// ---- resolve totally-split primes ----
printf "\nResolving totally-split primes with lambda...\n";
full := [];   // <Np, p, rid, a_p mod 11>
nclean := 0; nsplit := 0; nfail := 0;
for row in frobdata do
    Np := row[1]; p := row[2]; rid := row[3];
    if p in {2,11} or p eq 66179 then continue; end if;
    hbar := HbarFromCoeffs(row[4]);
    typ, lin, quad := Classify(hbar);
    Npbar := F11 ! Np;
    if typ eq "clean" then
        ap := -Coefficient(quad,1);
        Append(~full, <Np, p, rid, Integers()!ap>);
        nclean +:= 1;
    elif typ eq "split" then
        P := PrimeFromTag(p, rid);
        lam := chi(lam_img, P);
        lamp := Npbar / lam;
        // remove {lam, lamp} from the 4-root multiset
        roots := lin;
        ok := true;
        for r in [lam, lamp] do
            idx := Index(roots, r);
            if idx eq 0 then ok := false; break; end if;
            Remove(~roots, idx);
        end for;
        if ok and #roots eq 2 then
            ap := roots[1] + roots[2];
            Append(~full, <Np, p, rid, Integers()!ap>);
            nsplit +:= 1;
        else
            nfail +:= 1;
        end if;
    end if;
end for;
printf "  clean: %o,  resolved split: %o,  failed: %o,  total a_p: %o\n",
    nclean, nsplit, nfail, #full;

// write full table
fn := "mod11_apdata_full.txt";
PrintFile(fn, "// Curve A: 2-dim mod-11 rep W, FULL a_p (clean + lambda-resolved splits)" : Overwrite := true);
PrintFile(fn, Sprintf("// lambda order %o; rows <N(p), p, rid, a_p mod 11>", ord));
PrintFile(fn, Sprintf("aplist_full := %o;", full));
printf "Wrote %o (%o entries).\n", fn, #full;

// ---- bonus: determine lambda's conductor support (ramification at 2? 11? 66179?) ----
printf "\nConductor support of lambda (does a consistent character survive on smaller moduli?):\n";
for spec in [<"pbad*11*p2^4", pbad*(11*OK)*p2^4>,
             <"pbad*11      ", pbad*(11*OK)>,
             <"pbad*p2^4    ", pbad*p2^4>,
             <"11*p2^4      ", (11*OK)*p2^4>,
             <"pbad         ", pbad>,
             <"11           ", 11*OK>] do
    g := FindLambda(spec[2]);
    printf "  modulus %o : %o consistent character(s)\n", spec[1], #g;
end for;
