// goal2_mod5_search.m -- search the BFS RM-by-sqrt3 family for the moduli point
// (a:b:c) in P^2(Q(sqrt2)) whose Jacobian is B_f (LMFDB 2.2.8.1-5000.1-j).
//
// Filter: at split primes l in {7,17,23,31,41,47}, B_f has point-count trace
// s_l = {-2,6,14,6,-8,0} (from goal2_mod5_Bf_out.txt), the SAME at both primes
// above l.  For a candidate C_{a,b,c}/Q(sqrt2) we reduce mod each prime (both
// roots of 2 mod l) and require trace l+1-#C(F_l) = s_l for BOTH roots.
//
// Run: magma goal2_mod5_search.m
SetColumns(0);

// ---- build the family sextic symbolically over Q[a,b,c] ----
K3<z3> := CyclotomicField(3);
R3<A3,B3,C3> := PolynomialRing(K3, 3);
Qf3 := func< X,Y,Z | X^2 + X*(Y-Z) + (Y-Z)^2 - 3*X*Z >;
t3s := Qf3(A3, z3*B3, C3) * Qf3(A3, z3^2*B3, C3);   // rational-valued
RQ<a,b,c> := PolynomialRing(Rationals(), 3);
mons := Monomials(t3s); cofs := Coefficients(t3s);
t3Q := &+[ (Rationals()!cofs[i]) * Monomial(RQ, Exponents(mons[i])) : i in [1..#mons] ];

RQx<aa,bb,cc,X> := PolynomialRing(Rationals(), 4);
phi := hom< RQ -> RQx | aa, bb, cc >;
t3x := phi(t3Q);
Qf := func< XX,YY,ZZ | XX^2 + XX*(YY-ZZ) + (YY-ZZ)^2 - 3*XX*ZZ >;
q4x := Qf(cc, bb, aa);
H1 := bb*(aa-cc)*q4x*X^2 - q4x*X - 1;
lam1 := 4*aa*cc*t3x;
quart := aa^4 - 4*aa^3*bb - 8*aa^3*cc + 3*aa^2*bb^2 + 18*aa^2*bb*cc + 18*aa^2*cc^2
         + 2*aa*bb^3 - 9*aa*bb^2*cc - 12*aa*bb*cc^2 - 8*aa*cc^3 - 2*bb^4 - bb^3*cc
         + 2*bb*cc^3 + cc^4;
cub := 2*aa^3 - 3*aa^2*bb - 9*aa^2*cc + 9*aa*bb*cc + 6*aa*cc^2 + bb^3 - cc^3;
con := aa^2 - 2*aa*bb - 4*aa*cc + bb^2 + 7*bb*cc + cc^2;
G1 := bb*(aa-cc)^2*(aa^2-4*aa*cc-bb^2+cc^2)*q4x^2*X^3
      - (aa-cc)*q4x*quart*X^2 - q4x*cub*X - aa*con;
SEX := G1^2 + lam1*H1^3;     // in RQx; degree <=6 in X

// extract x-coefficients as polynomials in (aa,bb,cc)
coeffX := [ Rationals() ! 0 : i in [1..7] ];  // placeholder types replaced below
CoefX := [* *];
for i in [0..6] do
    Append(~CoefX, Coefficient(SEX, X, i));    // element of RQx in aa,bb,cc
end for;

// evaluate the 7 coeffs at numeric (a,b,c) in a field, returns list in that field
function SexCoeffs(av, bv, cv)
    return [ Evaluate(CoefX[i], [av, bv, cv, Parent(av)!0]) : i in [1..7] ];
end function;

// ---- fast point count of y^2 = sum cf[i] x^{i-1} over GF(l) ----
function FastCount(cf, l)
    Fl := GF(l);
    cc2 := [ Fl!x : x in cf ];
    d := 6; while d ge 0 and cc2[d+1] eq 0 do d -:= 1; end while;
    if d lt 3 then return -1; end if;        // degenerate reduction, skip
    N := 0;
    for xx in Fl do
        val := &+[ cc2[i]*xx^(i-1) : i in [1..d+1] ];
        if val eq 0 then N +:= 1;
        elif IsSquare(val) then N +:= 2; end if;
    end for;
    if IsEven(d) then
        if IsSquare(cc2[d+1]) then N +:= 2; end if;
    else
        N +:= 1;
    end if;
    return N;
end function;

// ---- targets & split-prime data ----
prs := [7,17,23,31,41,47];
strace := [-2, 6, 14, 6, -8, 0];
roots := AssociativeArray();
for l in prs do
    roots[l] := [ Integers()!r[1] : r in Roots(PolynomialRing(GF(l))![-2,0,1]) ];
end for;

K<s2> := QuadraticField(2);

// reduce a K-element e = p + q*s2 mod l with s2 -> r
function RedK(e, l, r)
    sq := Eltseq(K!e);      // [p, q], rationals
    return GF(l)!(Integers()!(Numerator(sq[1])) ) / GF(l)!(Integers()!Denominator(sq[1]))
         + (GF(l)!(Integers()!Numerator(sq[2]))/GF(l)!Denominator(sq[2])) * (GF(l)!r);
end function;

// |trace| check at prime l (both primes above l) against |s| -- twist-invariant,
// and does NOT force the two conjugate primes to share the sign of the trace.
function TraceOK(cfK, l, s)
    for r in roots[l] do
        cfl := [ RedK(cfK[i], l, r) : i in [1..7] ];
        N := FastCount(cfl, l);
        if N eq -1 then return false; end if;
        if AbsoluteValue(l + 1 - N) ne AbsoluteValue(s) then return false; end if;
    end for;
    return true;
end function;

// ---- search ----
B := 6;
vals := [ p + q*s2 : p in [-B..B], q in [-B..B] ];
LOGF := "goal2_mod5_search_out.txt";
PrintFile(LOGF, Sprintf("# BFS moduli search, height bound B=%o over Z[sqrt2]", B) : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;
LOG(Sprintf("candidates per coord: %o ; scanning ...", #vals));

cnt := 0; hits := 0; t0 := Cputime();
for av in vals do
    if av eq 0 then continue; end if;
    for cv in vals do
        if cv eq 0 or cv eq av then continue; end if;
        for bv in vals do
            if bv eq 0 then continue; end if;
            cnt +:= 1;
            cfK := SexCoeffs(av, bv, cv);
            // filter cheaply on l=7 first, then 17, 23, then the rest
            if not TraceOK(cfK, 7, strace[1]) then continue; end if;
            if not TraceOK(cfK, 17, strace[2]) then continue; end if;
            if not TraceOK(cfK, 23, strace[3]) then continue; end if;
            ok := true;
            for j in [4..6] do
                if not TraceOK(cfK, prs[j], strace[j]) then ok := false; break; end if;
            end for;
            if ok then
                hits +:= 1;
                LOG(Sprintf("HIT (a,b,c) = (%o, %o, %o)", av, bv, cv));
            end if;
        end for;
    end for;
end for;
LOG(Sprintf("scanned %o triples in %.1o s; %o hits", cnt, Cputime(t0), hits));
exit;
