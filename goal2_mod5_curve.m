// goal2_mod5_curve.m -- build the genus-2 curve for B_f from the reconstructed
// H_12 moduli point (e,f) = (-431/2341, 808/643) (rational!), and verify it against
// B_f's L-poly fingerprint over Q(sqrt2) up to quadratic twist.
SetColumns(0);
function IGtuple(e,f)
    A1:=f-1; A:=-(f^4+15*e*f+9*e)/3; B1:=(2*f^3-2*f^2-3*f+3*e+3)/3;
    B:=(2*f^6-63*e*f^3-81*e*f^2-54*e^2)/27; B2:=e^3;
    return [ -24*B1/A1, -12*A, 96*(A/A1)*B1-36*B, -4*A1*B2 ];
end function;

e0 := -431/2341; f0 := 808/643;
QQ := Rationals();
ic := IGtuple(e0,f0);
printf "Igusa-Clebsch invariants (rational): %o\n", ic;
C := HyperellipticCurveFromIgusaClebsch(ic);
print "=== Genus-2 curve C with RM by sqrt3 (moduli point of B_f) ===";
printf "base field of Mestre model: %o\n", BaseRing(C);
print C;

// Verify by reducing the (rational) Igusa-Clebsch invariants mod each split prime
// of Q(sqrt2) and matching B_f's L-poly up to quadratic twist.  (sqrt2 in F_l via
// each root; both primes above l checked.)
print "\n=== L-poly verification vs B_f fingerprint (up to twist) ===";
RTq<Tq> := PolynomialRing(QQ);
Ltgt := func< l,s,n | l^2*Tq^4 - l*s*Tq^3 + (2*l+n)*Tq^2 - s*Tq + 1 >;
tg := [ <7,-2,-2>,<17,6,-3>,<23,14,46>,<31,6,6>,<41,-8,-11>,<47,0,-48>,<71,4,-8>,
        <73,-8,13>,<79,2,-2>,<89,-6,-3>,<97,-4,-8>,<103,16,16>,<113,32,253>,
        <127,18,-66>,<137,-26,121>,<151,14,-98> ];
nok := 0; nbad := 0;
for t in tg do
    l := t[1]; PF<TT> := PolynomialRing(GF(l));
    Lt := PF!Ltgt(l,t[2],t[3]); Lw := PF!Ltgt(l,-t[2],t[3]);
    ok := true;
    for r in [Integers()!y[1] : y in Roots(PF![-2,0,1])] do   // sqrt2 -> r
        // reduce the rational Igusa-Clebsch invariants mod l (e,f are rational so
        // ic is rational; r only matters if ic had sqrt2 -- it doesn't -- but we
        // still check both primes for completeness)
        icl := [ GF(l)!x : x in ic ];
        if icl[4] eq 0 then ok := false; break; end if;
        try Cr := HyperellipticCurveFromIgusaClebsch(icl); catch err ok := false; break; end try;
        if LPolynomial(Cr) notin {Lt,Lw} then ok := false; break; end if;
    end for;
    if ok then nok +:= 1; else nbad +:= 1; printf "  MISMATCH at l=%o\n", l; end if;
end for;
printf "\nsplit primes matching B_f L-poly (up to twist): %o / %o\n", nok, nok+nbad;
if nbad eq 0 then
    print ">>> C reproduces B_f's fingerprint at ALL tested split primes.";
    print "    So Jac(C) is the RM-by-sqrt3 abelian surface of B_f (up to twist),";
    print "    with rational Igusa-Clebsch invariants (field of moduli Q).";
end if;
exit;
