// goal2_mod5_search2.m -- targeted search: since cond(B_f)=2^9*5^5, the model
// discriminant 2^12 a^3 b^9 c^14 (a-c)^12 q3^3 t3^3 q4^4 t4^4 must be supported at
// {2,5}, so a,b,c are {2,5}-S-units of Q(sqrt2).  Normalize c=1; scan a,b over a
// small S-unit set; filter by |trace| at split primes against B_f's fingerprint.
SetColumns(0);

// ---- family sextic symbolic over Q[a,b,c] ----
K3<z3> := CyclotomicField(3);
R3<A3,B3,C3> := PolynomialRing(K3, 3);
Qf3 := func< X,Y,Z | X^2 + X*(Y-Z) + (Y-Z)^2 - 3*X*Z >;
t3s := Qf3(A3, z3*B3, C3) * Qf3(A3, z3^2*B3, C3);
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
quart := aa^4-4*aa^3*bb-8*aa^3*cc+3*aa^2*bb^2+18*aa^2*bb*cc+18*aa^2*cc^2+2*aa*bb^3-9*aa*bb^2*cc-12*aa*bb*cc^2-8*aa*cc^3-2*bb^4-bb^3*cc+2*bb*cc^3+cc^4;
cub := 2*aa^3-3*aa^2*bb-9*aa^2*cc+9*aa*bb*cc+6*aa*cc^2+bb^3-cc^3;
con := aa^2-2*aa*bb-4*aa*cc+bb^2+7*bb*cc+cc^2;
G1 := bb*(aa-cc)^2*(aa^2-4*aa*cc-bb^2+cc^2)*q4x^2*X^3 - (aa-cc)*q4x*quart*X^2 - q4x*cub*X - aa*con;
SEX := G1^2 + lam1*H1^3;
CoefX := [* Coefficient(SEX, X, i) : i in [0..6] *];
function SexCoeffs(av,bv,cv) return [ Evaluate(CoefX[i], [av,bv,cv,Parent(av)!0]) : i in [1..7] ]; end function;

function FastCount(cf, l)
    Fl := GF(l); cc2 := [Fl!x : x in cf];
    d := 6; while d ge 0 and cc2[d+1] eq 0 do d -:= 1; end while;
    if d lt 3 then return -1; end if;
    N := 0;
    for xx in Fl do
        val := &+[cc2[i]*xx^(i-1) : i in [1..d+1]];
        if val eq 0 then N +:= 1; elif IsSquare(val) then N +:= 2; end if;
    end for;
    if IsEven(d) then if IsSquare(cc2[d+1]) then N +:= 2; end if; else N +:= 1; end if;
    return N;
end function;

prs := [7,17,23,31,41,47]; strace := [-2,6,14,6,-8,0];
roots := AssociativeArray();
for l in prs do roots[l] := [Integers()!r[1] : r in Roots(PolynomialRing(GF(l))![-2,0,1])]; end for;
K<s2> := QuadraticField(2);
function RedK(e, l, r)
    sq := Eltseq(K!e); Fl := GF(l);
    return (Fl!Numerator(sq[1]))/(Fl!Denominator(sq[1])) + ((Fl!Numerator(sq[2]))/(Fl!Denominator(sq[2])))*(Fl!r);
end function;
function TraceOK(cfK, l, s)
    for r in roots[l] do
        cfl := [RedK(cfK[i], l, r) : i in [1..7]];
        N := FastCount(cfl, l); if N eq -1 then return false; end if;
        if AbsoluteValue(l+1-N) ne AbsoluteValue(s) then return false; end if;
    end for;
    return true;
end function;

// ---- S-unit set for S={2,5} in Q(sqrt2): +-(1+s2)^i * s2^j * 5^k ----
u := 1 + s2;
suset := {K| };
for si in [1,-1] do for i in [-4..4] do for j in [-3..4] do for k in [-1..1] do
    Include(~suset, si*u^i*s2^j*5^k);
end for; end for; end for; end for;
su := Sort(SetToSequence(suset), func<x,y | AbsoluteValue(Norm(x)) le AbsoluteValue(Norm(y)) select -1 else 1>);
LOGF := "goal2_mod5_search2_out.txt";
PrintFile(LOGF, Sprintf("# S-unit search (c=1), #S-units=%o", #su) : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;
LOG(Sprintf("S-unit set size %o ; scanning a,b (c=1) ...", #su));

cv := K!1; cnt := 0; hits := 0; t0 := Cputime();
for av in su do
    if av eq cv then continue; end if;
    for bv in su do
        cnt +:= 1;
        cfK := SexCoeffs(av,bv,cv);
        if not TraceOK(cfK,7,strace[1]) then continue; end if;
        if not TraceOK(cfK,17,strace[2]) then continue; end if;
        if not TraceOK(cfK,23,strace[3]) then continue; end if;
        ok := true;
        for j in [4..6] do if not TraceOK(cfK,prs[j],strace[j]) then ok:=false; break; end if; end for;
        if ok then hits +:= 1; LOG(Sprintf("HIT (a,b,c)=(%o, %o, 1)", av, bv)); end if;
    end for;
end for;
LOG(Sprintf("scanned %o pairs in %.1o s; %o hits", cnt, Cputime(t0), hits));
exit;
