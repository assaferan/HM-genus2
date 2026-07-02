// refine_structure.m -- settle the structure of A[11] (2+1+1 vs 2+2) from raw Frobenius data.
//
// For each good prime, h_P = x^4 + c1 x^3 + c2 x^2 + c3 x + c4 (mod 11), c4 = N(P)^2.
// We compute factor types, test the "never both 2-dim pieces irreducible" condition, and
// run a careful search for a global character lambda (the 11-isogeny line) of order | 10.
//
// Run: magma refine_structure.m

SetColumns(0);
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11); R11<x> := PolynomialRing(F11);
load "mod11_frobdata.txt";   // frobdata := [ <N, p, rid, [c0..c4]>, ... ]

LOGF := "refine_structure_out.txt";
PrintFile(LOGF, "# refine_structure" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;

// ---- factor types + structural checks ----
types := AssociativeArray();
nbothirred := 0; ntot := 0;
n112_detok := 0; n112 := 0;
lin_data := [];   // for [1,1,2] primes: <ideal, {r1,r2}, N>  (linear-root pair, det)
for row in frobdata do
    N := row[1]; p := row[2]; rid := row[3]; c := row[4];
    Nm := F11!N;
    h := x^4 + (F11!c[2])*x^3 + (F11!c[3])*x^2 + (F11!c[4])*x + (F11!c[5]);
    fac := Factorization(h);
    degs := Sort(&cat[ [Degree(f[1]) : i in [1..f[2]]] : f in fac ]);
    key := Sprintf("%o", degs);
    if IsDefined(types, key) then types[key] +:= 1; else types[key] := 1; end if;
    ntot +:= 1;

    // "both 2-dim pieces irreducible" = h is a product of two DISTINCT irreducible quadratics
    quads := [f[1] : f in fac | Degree(f[1]) eq 2];
    nq := &+[Integers()| f[2] : f in fac | Degree(f[1]) eq 2];
    if nq eq 2 and #[f : f in fac | Degree(f[1]) eq 1] eq 0 then nbothirred +:= 1; end if;

    // [1,1,2]: two linear + one irreducible quadratic
    lins := &cat[ [Roots(f[1])[1][1] : i in [1..f[2]]] : f in fac | Degree(f[1]) eq 1 ];
    if #lins eq 2 and nq eq 1 then
        n112 +:= 1;
        r1 := lins[1]; r2 := lins[2];
        q := [f[1] : f in fac | Degree(f[1]) eq 2][1];   // monic x^2 - a x + b
        b := Coefficient(q,0);                            // det of the quadratic piece
        // 2+1+1 clean iff linear pair has det N and quadratic has det N
        if (r1*r2 eq Nm) and (b eq Nm) then n112_detok +:= 1; end if;
        Append(~lin_data, <PrimeFromTag(p,rid), {r1, r2}, Nm>);
    end if;
end for;
LOG("factor types over " * Sprintf("%o", ntot) * " good primes:");
for key in Sort([k : k in Keys(types)]) do LOG(Sprintf("  %o : %o", key, types[key])); end for;
LOG(Sprintf("both-pieces-irreducible (type [2,2]): %o   <-- 0 expected if 2+1+1 or correlated", nbothirred));
LOG(Sprintf("[1,1,2] primes: %o; of these, det-clean (lin prod = N AND quad const = N): %o",
    n112, n112_detok));

// ---- search for a global character lambda (the 11-isogeny line), order | 10 ----
// F_11^* is cyclic of order 10, generator g0 = 2. dlog base 2.
g0 := F11!2;
dlog := AssociativeArray();
v := F11!1; for i in [0..9] do dlog[v] := i; v := v*g0; end for;
function DL(a) return dlog[a]; end function;   // F_11^* -> Z/10

P66179 := [Pp[1] : Pp in Factorization(66179*OK) | Norm(Pp[1]) eq 66179][1];
p2 := Factorization(2*OK)[1][1];

LOG("lambda search (char of order | 10, lin root at every [1,1,2] prime):");
best := 0; bestdesc := "";
moduli := [
   <P66179, "P66179">, <P66179*p2^2, "P66179*p2^2">, <P66179*p2^4, "P66179*p2^4">,
   <P66179*ideal<OK|11>, "P66179*11">, <P66179*p2^4*ideal<OK|11>, "P66179*p2^4*11">
];
for mm in moduli do
  for inf in [[1,2],[Integers()|]] do
    if #inf eq 0 then Rc, mp := RayClassGroup(mm[1]); else Rc, mp := RayClassGroup(mm[1], inf); end if;
    // enumerate all homomorphisms Rc -> Z/10 via the dual; Z/10 ~ F11^* additively
    Z10 := AbelianGroup([10]);
    Homs := Homomorphisms(Rc, Z10);    // all characters of order dividing 10
    for phi in Homs do
        good := 0;
        for ld in lin_data do
            gp := ld[1] @@ mp;
            val := Eltseq(phi(gp))[1] mod 10;     // in Z/10 = dlog of lambda(P)
            // lambda(P) must be one of the two linear roots
            if val eq DL(SetToSequence(ld[2])[1]) or
               (#ld[2] eq 2 and val eq DL(SetToSequence(ld[2])[2])) or
               (#ld[2] eq 1 and val eq DL(SetToSequence(ld[2])[1])) then
                good +:= 1;
            end if;
        end for;
        if good gt best then best := good; bestdesc := Sprintf("modulus %o, oo=%o", mm[2], inf); end if;
    end for;
  end for;
end for;
LOG(Sprintf("  best global-character fit: %o / %o  (%o)", best, #lin_data, bestdesc));
if best eq #lin_data then
    LOG("  => a rational line EXISTS: structure is 2+1+1, target = W (mod11_apdata.txt).");
else
    LOG("  => NO rational line: structure is NOT 2+1+1; the never-[2,2] correlation has another source.");
end if;
LOG("DONE.");
exit;
