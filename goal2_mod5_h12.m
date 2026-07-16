// goal2_mod5_h12.m -- per-prime diagnostic for the Elkies-Kumar Humbert surface
// H_12 (RM by sqrt3, NO level structure).  Igusa-Clebsch [I2,I4,I6,I10](e,f) from
// arXiv:1209.3527 ancillary igusa12.txt:
//   A1=f-1; A=-(f^4+15ef+9e)/3; B1=(2f^3-2f^2-3f+3e+3)/3;
//   B=(2f^6-63ef^3-81ef^2-54e^2)/27; B2=e^3;
//   [I2,I4,I6,I10] = [-24 B1/A1, -12 A, 96 (A/A1) B1 - 36 B, -4 A1 B2].
// For each split prime l scan (e,f) in F_l^2, build the genus-2 curve from these
// invariants, and test whether its L-poly = B_f's target L-poly (up to quad twist).
// BFS gave 0 matches at l=7,23,31; if H_12 matches at ALL of 7,17,23,31 it is the
// right family and we can then CRT-reconstruct the moduli point over Q(sqrt2).
SetColumns(0);

function IGtuple(e, f)
    A1 := f-1;
    A  := -(f^4+15*e*f+9*e)/3;
    B1 := (2*f^3-2*f^2-3*f+3*e+3)/3;
    B  := (2*f^6-63*e*f^3-81*e*f^2-54*e^2)/27;
    B2 := e^3;
    return [ -24*B1/A1, -12*A, 96*(A/A1)*B1-36*B, -4*A1*B2 ], A1;
end function;

badlocusf := func< e,f | e*f*(f^2-1)*(e-f^2)*(8*f^4-9*e*f^2-3*e)*(f^6-f^4-18*e*f^2+27*e^2+16*e) >;

RTq<Tq> := PolynomialRing(Rationals());
function Ltarget(l,s,n) return l^2*Tq^4 - l*s*Tq^3 + (2*l+n)*Tq^2 - s*Tq + 1; end function;
target := AssociativeArray();
target[7] := <-2,-2>; target[17] := <6,-3>; target[23] := <14,46>; target[31] := <6,6>;

for l in [7,17,23,31] do
    Fl := GF(l); PF<TT> := PolynomialRing(Fl);
    sn := target[l];
    Lt := PF ! Ltarget(l, sn[1], sn[2]);
    Lw := PF ! Ltarget(l, -sn[1], sn[2]);
    npts := 0; igset := {};
    for ev in Fl do for fv in Fl do
        if badlocusf(ev,fv) eq 0 then continue; end if;
        ic, A1 := IGtuple(ev, fv);
        if ic[4] eq 0 then continue; end if;         // I10=0 -> not genus 2
        try
            C := HyperellipticCurveFromIgusaClebsch(ic);
        catch err
            continue;
        end try;
        L := LPolynomial(C);
        if L eq Lt or L eq Lw then
            npts +:= 1;
            Include(~igset, G2Invariants(C));
        end if;
    end for; end for;
    printf "l=%o : %o matching (e,f) points, %o distinct Igusa tuples\n", l, npts, #igset;
end for;
exit;
