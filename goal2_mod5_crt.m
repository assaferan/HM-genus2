// goal2_mod5_crt.m -- reconstruct the H_12 moduli point (e,f) in Q(sqrt2)^2 of B_f
// by CRT over split primes.  At each split prime l the matching (e,f) in F_l^2 is a
// small set; the true (e0,f0) reduces at the two primes above l (sqrt2->r1,r2) into
// it.  From a conjugate pair we get (e0,f0) mod l in Q(sqrt2); CRT + rational
// reconstruction over primes 7,23,17,31,41 (modulus ~3.5e6, height up to ~1300).
// Verify survivors at fresh split primes 71,73; then build the curve over Q(sqrt2).
SetColumns(0);

function IGtuple(e, f)
    A1 := f-1; A := -(f^4+15*e*f+9*e)/3; B1 := (2*f^3-2*f^2-3*f+3*e+3)/3;
    B := (2*f^6-63*e*f^3-81*e*f^2-54*e^2)/27; B2 := e^3;
    return [ -24*B1/A1, -12*A, 96*(A/A1)*B1-36*B, -4*A1*B2 ];
end function;
badlocusf := func< e,f | e*f*(f^2-1)*(e-f^2)*(8*f^4-9*e*f^2-3*e)*(f^6-f^4-18*e*f^2+27*e^2+16*e) >;
RTq<Tq> := PolynomialRing(Rationals());
function Ltgt(l,s,n) return l^2*Tq^4 - l*s*Tq^3 + (2*l+n)*Tq^2 - s*Tq + 1; end function;
target := AssociativeArray();
for t in [ <7,-2,-2>,<17,6,-3>,<23,14,46>,<31,6,6>,<41,-8,-11>,<71,4,-8>,<73,-8,13> ] do
    target[t[1]] := <t[2],t[3]>;
end for;

function MatchEF(l)
    Fl := GF(l); PF<TT> := PolynomialRing(Fl); sn := target[l];
    Lt := PF ! Ltgt(l,sn[1],sn[2]); Lw := PF ! Ltgt(l,-sn[1],sn[2]);
    M := [];
    for ev in Fl do for fv in Fl do
        if badlocusf(ev,fv) eq 0 then continue; end if;
        ic := IGtuple(ev,fv); if ic[4] eq 0 then continue; end if;
        try C := HyperellipticCurveFromIgusaClebsch(ic); catch e continue; end try;
        L := LPolynomial(C);
        if L eq Lt or L eq Lw then Append(~M, <ev,fv>); end if;
    end for; end for;
    return M;
end function;

LOGF := "goal2_mod5_crt_out.txt";
PrintFile(LOGF, "# CRT reconstruction of B_f moduli point (e,f) on H_12" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

primes := [7,23,31,41,73];
Mset := AssociativeArray(); rts := AssociativeArray();
for l in primes do
    Mset[l] := MatchEF(l);
    rts[l] := [Integers()!x[1] : x in Roots(PolynomialRing(GF(l))![-2,0,1])];
    LOG(Sprintf("l=%o: |M|=%o", l, #Mset[l]));
end for;

// per-prime candidate residues (xe,ye,xf,yf) in Z from an ordered pair (p1,p2):
// p1 = (e0,f0) mod (sqrt2->r1), p2 = mod (sqrt2->r2).
function CandsAt(l)
    Fl := GF(l); r1 := Fl!rts[l][1]; r2 := Fl!rts[l][2]; M := Mset[l];
    out := [];
    for p1 in M do for p2 in M do
        ye := (Fl!p1[1]-Fl!p2[1])/(r1-r2); xe := Fl!p1[1]-ye*r1;
        yf := (Fl!p1[2]-Fl!p2[2])/(r1-r2); xf := Fl!p1[2]-yf*r1;
        Append(~out, [Integers()!xe,Integers()!ye,Integers()!xf,Integers()!yf]);
    end for; end for;
    return out;
end function;

CA := [ CandsAt(l) : l in primes ];
Mprod := &*primes; ZM := Integers(Mprod);
LOG(Sprintf("combining %o candidate-tuples (modulus %o) ...", &*[#c : c in CA], Mprod));

sols := {}; np := #primes;
for c1 in CA[1] do for c2 in CA[2] do for c3 in CA[3] do for c4 in CA[4] do for c5 in CA[5] do
    comps := [c1,c2,c3,c4,c5];
    // reconstruct xe first (prune)
    vxe := CRT([comps[i][1] : i in [1..np]], primes);
    b1, qxe := RationalReconstruction(ZM!vxe); if not b1 then continue; end if;
    vye := CRT([comps[i][2] : i in [1..np]], primes);
    b2, qye := RationalReconstruction(ZM!vye); if not b2 then continue; end if;
    vxf := CRT([comps[i][3] : i in [1..np]], primes);
    b3, qxf := RationalReconstruction(ZM!vxf); if not b3 then continue; end if;
    vyf := CRT([comps[i][4] : i in [1..np]], primes);
    b4, qyf := RationalReconstruction(ZM!vyf); if not b4 then continue; end if;
    Include(~sols, <qxe,qye,qxf,qyf>);
end for; end for; end for; end for; end for;
LOG(Sprintf("reconstructed %o candidate (e,f); verifying at fresh primes 71,73 ...", #sols));

K<s2> := QuadraticField(2);
// precompute matching sets at fresh primes 71,73 for O(1) membership filtering
fprimes := [17,71];
Mfset := AssociativeArray(); frts := AssociativeArray();
for l in fprimes do
    Mfset[l] := { <x[1],x[2]> : x in MatchEF(l) };
    frts[l] := [Integers()!y[1] : y in Roots(PolynomialRing(GF(l))![-2,0,1])];
    LOG(Sprintf("  |M_%o|=%o", l, #Mfset[l]));
end for;
function RedKlSafe(el, l, r)
    sq := Eltseq(K!el); Fl := GF(l);
    d1 := Denominator(sq[1]); d2 := Denominator(sq[2]);
    if (d1 mod l eq 0) or (d2 mod l eq 0) then return false, Fl!0; end if;
    return true, (Fl!Numerator(sq[1]))/(Fl!d1) + ((Fl!Numerator(sq[2]))/(Fl!d2))*(Fl!r);
end function;
function PassPrime(e0, f0, l)
    for r in frts[l] do
        o1,eb := RedKlSafe(e0,l,r); o2,fb := RedKlSafe(f0,l,r);
        if (not o1) or (not o2) then return false; end if;
        if <eb,fb> notin Mfset[l] then return false; end if;
    end for;
    return true;
end function;

verified := [];
for q in sols do
    e0 := q[1]+q[2]*s2; f0 := q[3]+q[4]*s2;
    if PassPrime(e0,f0,17) and PassPrime(e0,f0,71) then
        Append(~verified, <e0,f0>);
        LOG(Sprintf("VERIFIED (e,f) = (%o, %o)", e0, f0));
    end if;
end for;

LOG(Sprintf("\n=== %o verified survivors; building curve over Q(sqrt2) ===", #verified));
seen := {};
for ef in verified do
    ic := IGtuple(ef[1],ef[2]);
    try C := HyperellipticCurveFromIgusaClebsch(ic); catch e continue; end try;
    try C := SimplifiedModel(C); catch e ; end try;
    ig := IgusaClebschInvariants(C);
    if ig in seen then continue; end if; Include(~seen, ig);
    LOG(Sprintf("(e,f) = (%o, %o)", ef[1], ef[2]));
    LOG(Sprintf("  C : %o", C));
    LOG(Sprintf("  hyperelliptic poly f: %o", HyperellipticPolynomials(C)));
end for;
exit;
