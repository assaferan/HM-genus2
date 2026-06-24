// identify_twist.m -- identify the quadratic twist eps from the clean (p,rid,eps) data
// saved by extract_eigensystems.m (eps_data.txt). Cheap (no Hilbert modular forms).
//
// Run: magma identify_twist.m

SetColumns(0);
K<w> := QuadraticField(5); OK := Integers(K);
load "eps_data.txt";   // epsraw := [ <p, rid, eps>, ... ]
function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;
pts := [ <PrimeFromTag(e[1],e[2]), e[3], e[1], e[2]> : e in epsraw ];
printf "%o data points.\n\n", #pts;
p2 := Factorization(2*OK)[1][1];

// best order-2 character of RayClassGroup(m,inf): returns <best count, failing pts>
function BestChar(m, inf)
    Rcl, mp := RayClassGroup(m, inf); ng := Ngens(Rcl);
    ords := [ Order(Rcl.i) : i in [1..ng] ];
    choices := [ (IsEven(ords[i]) select [1,-1] else [1]) : i in [1..ng] ];
    best := -1; bestfail := [];
    for tup in CartesianProduct(choices) do
        img := [ tup[i] : i in [1..ng] ]; c := 0; fail := [];
        for pe in pts do
            es := Eltseq(pe[1] @@ mp);
            if &*[ img[i]^es[i] : i in [1..ng] ] eq pe[2] then c +:= 1;
            else Append(~fail, <pe[3],pe[4]>); end if;
        end for;
        if c gt best then best := c; bestfail := fail; end if;
    end for;
    return best, bestfail;
end function;

for ms in [<"(2)^4 (5) (66179)",       p2^4*(5*OK)*(66179*OK)>,
           <"(2)^4 (5) (11) (66179)",  p2^4*(5*OK)*(11*OK)*(66179*OK)>,
           <"(2)^8 (5) (66179)",       p2^8*(5*OK)*(66179*OK)>,
           <"(2)^4 (5) (66179)^2",     p2^4*(5*OK)*(66179*OK)^2>] do
    for inf in [[1,2],[Integers()|]] do
        b, fl := BestChar(ms[2], inf);
        printf "modulus %-24o inf=%-7o: best %o/%o; fails at %o\n",
            ms[1], inf, b, #pts, fl;
    end for;
end for;
