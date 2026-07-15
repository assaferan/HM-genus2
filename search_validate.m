// validate counting machinery vs Magma on the sample C_{1,2,-1}.
SetColumns(0);
K3<z3> := CyclotomicField(3);
Qf := func< X,Y,Z | X^2 + X*(Y-Z) + (Y-Z)^2 - 3*X*Z >;
function t3val(a,b,c) return Rationals()!(Qf(a,z3*b,c)*Qf(a,z3^2*b,c)); end function;
function BFSsextic(a,b,c,Px)
    x := Px.1; q4 := Qf(c,b,a); t3 := t3val(a,b,c);
    H1 := b*(a-c)*q4*x^2 - q4*x - 1; lam1 := 4*a*c*t3;
    quart := a^4-4*a^3*b-8*a^3*c+3*a^2*b^2+18*a^2*b*c+18*a^2*c^2+2*a*b^3-9*a*b^2*c-12*a*b*c^2-8*a*c^3-2*b^4-b^3*c+2*b*c^3+c^4;
    cub := 2*a^3-3*a^2*b-9*a^2*c+9*a*b*c+6*a*c^2+b^3-c^3;
    con := a^2-2*a*b-4*a*c+b^2+7*b*c+c^2;
    G1 := b*(a-c)^2*(a^2-4*a*c-b^2+c^2)*q4^2*x^3 - (a-c)*q4*quart*x^2 - q4*cub*x - a*con;
    return G1^2 + lam1*H1^3;
end function;
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

QQ := Rationals(); Px<x> := PolynomialRing(QQ);
fam := BFSsextic(1,2,-1,Px);
cf := [ Coefficient(fam, i) : i in [0..6] ];
print "sample family sextic:", fam;
C := HyperellipticCurve(fam);
for l in [7,17,23,31,41,47] do
    magcnt := #ChangeRing(C, GF(l));
    fastcnt := FastCount([GF(l)!c : c in cf], l);
    printf "l=%o : Magma #C=%o  Fast #C=%o  trace=%o  %o\n",
        l, magcnt, fastcnt, l+1-magcnt, magcnt eq fastcnt select "OK" else "MISMATCH";
end for;
exit;
