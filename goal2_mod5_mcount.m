// count |M_l| (matching (e,f) up to twist) for split primes, to choose CRT primes.
SetColumns(0);
function IGtuple(e,f)
    A1:=f-1; A:=-(f^4+15*e*f+9*e)/3; B1:=(2*f^3-2*f^2-3*f+3*e+3)/3;
    B:=(2*f^6-63*e*f^3-81*e*f^2-54*e^2)/27; B2:=e^3;
    return [ -24*B1/A1, -12*A, 96*(A/A1)*B1-36*B, -4*A1*B2 ];
end function;
badl := func< e,f | e*f*(f^2-1)*(e-f^2)*(8*f^4-9*e*f^2-3*e)*(f^6-f^4-18*e*f^2+27*e^2+16*e) >;
RTq<Tq> := PolynomialRing(Rationals());
Lt := func< l,s,n | l^2*Tq^4 - l*s*Tq^3 + (2*l+n)*Tq^2 - s*Tq + 1 >;
// (s,n) targets
tgt := [ <7,-2,-2>,<17,6,-3>,<23,14,46>,<31,6,6>,<41,-8,-11>,<47,0,-48>,
         <71,4,-8>,<73,-8,13>,<79,2,-2>,<89,-6,-3>,<97,-4,-8> ];
for t in tgt do
    l:=t[1]; s:=t[2]; n:=t[3];
    Fl:=GF(l); PF<TT>:=PolynomialRing(Fl);
    L1:=PF!Lt(l,s,n); L2:=PF!Lt(l,-s,n);
    cnt:=0;
    for ev in Fl do for fv in Fl do
        if badl(ev,fv) eq 0 then continue; end if;
        ic:=IGtuple(ev,fv); if ic[4] eq 0 then continue; end if;
        try C:=HyperellipticCurveFromIgusaClebsch(ic); catch e continue; end try;
        L:=LPolynomial(C);
        if L eq L1 or L eq L2 then cnt+:=1; end if;
    end for; end for;
    printf "l=%o |M|=%o\n", l, cnt;
end for;
exit;
