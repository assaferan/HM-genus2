// identify_twist2.m -- find eps directly in the Hecke character group, exclude degenerate
// (double-root) primes, and read off its conductor and quadratic field K(sqrt d).
//
// Run: magma identify_twist2.m

SetColumns(0);
K<w> := QuadraticField(5); OK := Integers(K);
load "eps_data.txt";    // epsraw := [ <p, rid, eps>, ... ]
load "mod11_trnm.txt";  // trnm
F11 := GF(11);

function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;

// degenerate primes: q has a double root (Tr^2 - 4 Nm = 0 mod 11) -> eps ill-defined
degen := {};
for r in trnm do
    if (F11!r[4])^2 - 4*(F11!r[5]) eq 0 then Include(~degen, <r[2],r[3]>); end if;
end for;
pts := [ <PrimeFromTag(e[1],e[2]), e[3], e[1], e[2]> : e in epsraw | <e[1],e[2]> notin degen ];
printf "%o data points, %o non-degenerate used.\n\n", #epsraw, #pts;

p2 := Factorization(2*OK)[1][1];
m := p2^4*(5*OK)*(66179*OK);
G := HeckeCharacterGroup(m, [1,2]);
printf "Hecke character group order %o; searching order-<=2 characters...\n", #G;

best := -1; bestchi := G!1;
for chi in Elements(G) do
    if Order(chi) gt 2 then continue; end if;
    if Order(chi) eq 1 then continue; end if;
    c := 0;
    for pe in pts do
        v := chi(pe[1]);              // +-1 for an order-2 character
        sign := (v eq 1) select 1 else -1;
        if sign eq pe[2] then c +:= 1; end if;
    end for;
    if c gt best then best := c; bestchi := chi; end if;
end for;
printf "best order-2 character matches %o / %o non-degenerate primes\n", best, #pts;

cond := Conductor(bestchi);
printf "  conductor of eps: %o  (norm %o)\n", Factorization(cond), Norm(cond);

// realize as K(sqrt d): the fixed field of ker(eps)
A := AbelianExtension(bestchi);
L := NumberField(A);
printf "  fixed field defining polynomial over K: %o\n", DefiningPolynomial(L);
// extract d from x^2 - d (after completing the square / simplifying)
OL := MaximalOrder(AbsoluteField(L));
printf "  K(sqrt d), discriminant data: disc(L/K) factorization = %o\n", Factorization(Discriminant(A));
