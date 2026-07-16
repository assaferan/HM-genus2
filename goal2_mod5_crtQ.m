// goal2_mod5_crtQ.m -- reconstruct B_f's H_12 moduli point assuming it is RATIONAL
// (field of moduli Q, forced by f's inner twist).  Then (e0,f0) mod l is a SINGLE
// point of M_l, so combos = prod|M_l| (not squared).  Use many primes -> big
// modulus, then verify at many fresh primes.
SetColumns(0);
function IGtuple(e,f)
    A1:=f-1; A:=-(f^4+15*e*f+9*e)/3; B1:=(2*f^3-2*f^2-3*f+3*e+3)/3;
    B:=(2*f^6-63*e*f^3-81*e*f^2-54*e^2)/27; B2:=e^3;
    return [ -24*B1/A1, -12*A, 96*(A/A1)*B1-36*B, -4*A1*B2 ];
end function;
badl := func< e,f | e*f*(f^2-1)*(e-f^2)*(8*f^4-9*e*f^2-3*e)*(f^6-f^4-18*e*f^2+27*e^2+16*e) >;
RTq<Tq> := PolynomialRing(Rationals());
Ltgt := func< l,s,n | l^2*Tq^4 - l*s*Tq^3 + (2*l+n)*Tq^2 - s*Tq + 1 >;
target := AssociativeArray();
for t in [ <7,-2,-2>,<17,6,-3>,<23,14,46>,<31,6,6>,<41,-8,-11>,<47,0,-48>,<71,4,-8>,
           <73,-8,13>,<79,2,-2>,<89,-6,-3>,<97,-4,-8>,<103,16,16>,<113,32,253>,
           <127,18,-66>,<137,-26,121>,<151,14,-98> ] do target[t[1]]:=<t[2],t[3]>; end for;

function MatchEF(l)
    Fl:=GF(l); PF<TT>:=PolynomialRing(Fl); sn:=target[l];
    Lt:=PF!Ltgt(l,sn[1],sn[2]); Lw:=PF!Ltgt(l,-sn[1],sn[2]); M:=[];
    for ev in Fl do for fv in Fl do
        if badl(ev,fv) eq 0 then continue; end if;
        ic:=IGtuple(ev,fv); if ic[4] eq 0 then continue; end if;
        try C:=HyperellipticCurveFromIgusaClebsch(ic); catch e continue; end try;
        if LPolynomial(C) in {Lt,Lw} then Append(~M,<ev,fv>); end if;
    end for; end for; return M;
end function;

LOGF:="goal2_mod5_crtQ_out.txt";
PrintFile(LOGF,"# rational CRT reconstruction of B_f moduli point":Overwrite:=true);
procedure LOG(s) printf "%o\n",s; PrintFile(LOGF,s); end procedure;

primes := [7,23,17,31,41,73,79];   // prod|M| ~ 62400*30, modulus ~2e10, H<1e5
Mset := AssociativeArray();
for l in primes do Mset[l]:=MatchEF(l); LOG(Sprintf("l=%o |M|=%o",l,#Mset[l])); end for;
Mprod := &*primes; ZM := Integers(Mprod);
LOG(Sprintf("combining %o tuples, modulus %o (height < %o) ...",
    &*[#Mset[l]:l in primes], Mprod, Isqrt(Mprod div 2)));

// enumerate one point per prime; CRT e-coords and f-coords; reconstruct
sols := {};
np := #primes;
Ms := [* Mset[l] : l in primes *];
sizes := [ #Ms[j] : j in [1..np] ];
total := &*sizes;
for cc in [0..total-1] do
    // decode mixed-radix index
    t := cc; sel := [];
    for j in [1..np] do Append(~sel, (t mod sizes[j])+1); t := t div sizes[j]; end for;
    eres := [ Integers()!(Ms[j][sel[j]][1]) : j in [1..np] ];
    fres := [ Integers()!(Ms[j][sel[j]][2]) : j in [1..np] ];
    be, qe := RationalReconstruction(ZM!CRT(eres, primes)); if not be then continue; end if;
    bf, qf := RationalReconstruction(ZM!CRT(fres, primes)); if not bf then continue; end if;
    Include(~sols, <qe,qf>);
end for;
LOG(Sprintf("reconstructed %o rational (e,f); fast membership filter at 47,79 ...", #sols));

// fast filter: precompute M_47, M_79 as sets; (e0 mod l, f0 mod l) must be in them
Mf := AssociativeArray();
for l in [47,89] do Mf[l] := { x : x in MatchEF(l) }; LOG(Sprintf("  |M_%o|=%o",l,#Mf[l])); end for;
function memOK(e0,f0,l)
    if (Denominator(e0) mod l eq 0) or (Denominator(f0) mod l eq 0) then return false; end if;
    return <GF(l)!e0, GF(l)!f0> in Mf[l];
end function;
cand2 := [ ef : ef in sols | memOK(ef[1],ef[2],47) and memOK(ef[1],ef[2],89) ];
LOG(Sprintf("%o pass 47 & 89; confirming at 89,97,103,113,127,137,151 ...", #cand2));

// confirm survivors at more fresh primes by building the curve mod l
verified := [];
for ef in cand2 do
    e0:=ef[1]; f0:=ef[2]; good:=true;
    for l in [97,103,113,127,137,151] do
        Fl:=GF(l); PF<TT>:=PolynomialRing(Fl); sn:=target[l];
        Lt:=PF!Ltgt(l,sn[1],sn[2]); Lw:=PF!Ltgt(l,-sn[1],sn[2]);
        if (Denominator(e0) mod l eq 0) or (Denominator(f0) mod l eq 0) then good:=false; break; end if;
        if badl(Fl!e0,Fl!f0) eq 0 then good:=false; break; end if;
        ic:=IGtuple(Fl!e0,Fl!f0); if ic[4] eq 0 then good:=false; break; end if;
        try C:=HyperellipticCurveFromIgusaClebsch([Fl|x:x in ic]); catch e good:=false; break; end try;
        if LPolynomial(C) notin {Lt,Lw} then good:=false; break; end if;
    end for;
    if good then Append(~verified,ef); LOG(Sprintf("VERIFIED (e,f)=(%o,%o)",e0,f0)); end if;
end for;
LOG(Sprintf("=== %o survivors passing ALL 16 fresh primes ===", #verified));
for ef in verified do LOG(Sprintf("  (e,f) = (%o, %o)\n  Igusa-Clebsch = %o", ef[1], ef[2], IGtuple(ef[1],ef[2]))); end for;
exit;
