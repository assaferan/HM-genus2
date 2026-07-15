// goal2_mod5_igusaset.m -- first stage of reduction-CRT recognition of B_f.
// At each small split prime l, scan (a:b:c) in P^2(F_l); keep those whose
// Jacobian L-poly equals B_f's target L-poly at l UP TO QUADRATIC TWIST (i.e.
// equals L_target or its s->-s twin); collect the distinct absolute-Igusa (G2)
// invariant tuples.  B_f's Igusa mod P is one of them (twist preserves Igusa).
// If the set J_l is small, CRT reconstruction of B_f's Igusa over Q(sqrt2) is easy.
SetColumns(0);

// symbolic family sextic over Q[a,b,c]
K3<z3> := CyclotomicField(3);
R3<A3,B3,C3> := PolynomialRing(K3, 3);
Qf3 := func< X,Y,Z | X^2 + X*(Y-Z) + (Y-Z)^2 - 3*X*Z >;
t3s := Qf3(A3, z3*B3, C3) * Qf3(A3, z3^2*B3, C3);
RQ<a,b,c> := PolynomialRing(Rationals(), 3);
mons := Monomials(t3s); cofs := Coefficients(t3s);
t3Q := &+[ (Rationals()!cofs[i]) * Monomial(RQ, Exponents(mons[i])) : i in [1..#mons] ];
RQx<aa,bb,cc,X> := PolynomialRing(Rationals(), 4);
phi := hom< RQ -> RQx | aa, bb, cc >; t3x := phi(t3Q);
Qf := func< XX,YY,ZZ | XX^2 + XX*(YY-ZZ) + (YY-ZZ)^2 - 3*XX*ZZ >;
q4x := Qf(cc, bb, aa);
H1 := bb*(aa-cc)*q4x*X^2 - q4x*X - 1; lam1 := 4*aa*cc*t3x;
quart := aa^4-4*aa^3*bb-8*aa^3*cc+3*aa^2*bb^2+18*aa^2*bb*cc+18*aa^2*cc^2+2*aa*bb^3-9*aa*bb^2*cc-12*aa*bb*cc^2-8*aa*cc^3-2*bb^4-bb^3*cc+2*bb*cc^3+cc^4;
cub := 2*aa^3-3*aa^2*bb-9*aa^2*cc+9*aa*bb*cc+6*aa*cc^2+bb^3-cc^3;
con := aa^2-2*aa*bb-4*aa*cc+bb^2+7*bb*cc+cc^2;
G1 := bb*(aa-cc)^2*(aa^2-4*aa*cc-bb^2+cc^2)*q4x^2*X^3 - (aa-cc)*q4x*quart*X^2 - q4x*cub*X - aa*con;
SEX := G1^2 + lam1*H1^3;
CoefX := [* Coefficient(SEX, X, i) : i in [0..6] *];

// target L-polys at split primes (from goal2_mod5_Bf_out.txt): (s, n)
target := AssociativeArray();
target[7]  := <-2, -2>;
target[17] := <6, -3>;
target[23] := <14, 46>;
target[31] := <6, 6>;

RT<T> := PolynomialRing(Rationals());
function Ltarget(l, s, n) return l^2*T^4 - l*s*T^3 + (2*l+n)*T^2 - s*T + 1; end function;

for l in [7,17,23,31] do
    Fl := GF(l); PF<TT> := PolynomialRing(Fl);
    sn := target[l]; s := sn[1]; n := sn[2];
    Lt := PF ! Ltarget(l, s, n);      // target L-poly mod l
    Lw := PF ! Ltarget(l, -s, n);     // its quadratic twist (s -> -s)
    igset := {};                       // distinct G2 invariant tuples
    npts := 0;
    // scan P^2(F_l): normalize first nonzero coord to 1
    reps := [];
    for x in Fl do Append(~reps, [Fl|1, x, 0]); end for;              // c=0 chart partial
    // full projective scan: (a:b:c), dedup by scaling
    seen := {};
    for av in Fl do for bv in Fl do for cv in Fl do
        if av eq 0 and bv eq 0 and cv eq 0 then continue; end if;
        // normalize
        v := [av,bv,cv];
        j := 1; while v[j] eq 0 do j +:= 1; end while;
        nv := [ x/v[j] : x in v ];
        if nv in seen then continue; end if; Include(~seen, nv);
        a0 := nv[1]; b0 := nv[2]; c0 := nv[3];
        if a0 eq 0 or c0 eq 0 or b0 eq 0 or a0 eq c0 then continue; end if;
        cf := [ Evaluate(CoefX[i], [a0,b0,c0,Fl!0]) : i in [1..7] ];
        d := 6; while d ge 0 and cf[d+1] eq 0 do d -:= 1; end while;
        if d lt 5 then continue; end if;
        fx := &+[ cf[i]*TT^(i-1) : i in [1..d+1] ];
        if Discriminant(fx) eq 0 then continue; end if;
        C := HyperellipticCurve(fx);
        L := LPolynomial(C);
        if L eq Lt or L eq Lw then
            npts +:= 1;
            Include(~igset, G2Invariants(C));
        end if;
    end for; end for; end for;
    printf "l=%o : %o matching points, %o DISTINCT Igusa tuples\n", l, npts, #igset;
    for ig in igset do printf "    %o\n", ig; end for;
end for;
exit;
