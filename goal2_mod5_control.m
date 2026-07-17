// goal2_mod5_control.m -- POSITIVE CONTROL for the H_12 CRT reconstruction pipeline.
// Take a KNOWN (e0,f0) in Q^2, build its H_12 curve, compute its L-poly fingerprint <s,n>
// at split primes, then run the SAME MatchEFfast + CRT reconstruction and check we recover
// (e0,f0).  Validates: (1) the true residue is in M_l; (2) CRT+RationalReconstruction works.
// Also reports the modulus/height-reach needed as a function of the control's height.
//
// Args: E, F (numerator forms as strings), e.g. E:="7/3" F:="-5/2".  Default (3,-2).
// Run: magma E:="7/3" F:="-5/2" goal2_mod5_control.m
SetColumns(0);
QQ := Rationals();
if assigned Es then e0 := eval("QQ!(" cat Es cat ")"); else e0 := QQ!3; end if;
if assigned Fs then f0 := eval("QQ!(" cat Fs cat ")"); else f0 := QQ!(-2); end if;

function IGtuple(e,f)
    A1:=f-1; A:=-(f^4+15*e*f+9*e)/3; B1:=(2*f^3-2*f^2-3*f+3*e+3)/3;
    B:=(2*f^6-63*e*f^3-81*e*f^2-54*e^2)/27; B2:=e^3;
    return [ -24*B1/A1, -12*A, 96*(A/A1)*B1-36*B, -4*A1*B2 ];
end function;
badl := func< e,f | e*f*(f^2-1)*(e-f^2)*(8*f^4-9*e*f^2-3*e)*(f^6-f^4-18*e*f^2+27*e^2+16*e) >;
function CountC(fp, Fq)
    N := 0; for x in Fq do v := Evaluate(fp, x);
        if v eq 0 then N +:= 1; elif IsSquare(v) then N +:= 2; end if; end for;
    d := Degree(fp);
    if IsEven(d) then if IsSquare(LeadingCoefficient(fp)) then N +:= 2; end if; else N +:= 1; end if;
    return N;
end function;
function MatchEFfast(l, s, n)
    Fl := GF(l); Fl2 := GF(l^2); M := [];
    for ev in Fl do for fv in Fl do
        if badl(ev,fv) eq 0 then continue; end if;
        ic := IGtuple(ev,fv); if ic[4] eq 0 then continue; end if;
        try C := HyperellipticCurveFromIgusaClebsch(ic); catch e continue; end try;
        fp := HyperellipticPolynomials(C); if Degree(fp) lt 5 then continue; end if;
        a1 := l + 1 - CountC(fp, Fl);
        if AbsoluteValue(a1) ne AbsoluteValue(s) then continue; end if;
        c2 := l^2 + 1 - CountC(ChangeRing(fp, Fl2), Fl2);
        e2 := (a1^2 - c2) div 2;
        if e2 - 2*l eq n then Append(~M, <Integers()!ev, Integers()!fv>); end if;
    end for; end for;
    return M;
end function;

LOGF := "goal2_mod5_control_out.txt";
PrintFile(LOGF, Sprintf("# control (e0,f0)=(%o,%o)  ht=%o", e0, f0, Max([Height(e0),Height(f0)])) : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;
LOG(Sprintf("control point (e0,f0) = (%o, %o), heights %o, %o", e0, f0, Height(e0), Height(f0)));

// fingerprint of the control curve at a prime l (reduce IGtuple mod l; get <s,n> from L-poly)
function ctrlSN(l)
    Fl := GF(l);
    if (Denominator(e0) mod l eq 0) or (Denominator(f0) mod l eq 0) then return false,0,0; end if;
    if badl(Fl!e0,Fl!f0) eq 0 then return false,0,0; end if;
    ic := IGtuple(Fl!e0,Fl!f0); if ic[4] eq 0 then return false,0,0; end if;
    try C := HyperellipticCurveFromIgusaClebsch([Fl|x:x in ic]); catch err return false,0,0; end try;
    fp := HyperellipticPolynomials(C); if Degree(fp) lt 5 then return false,0,0; end if;
    a1 := l+1 - CountC(fp, Fl);
    c2 := l^2+1 - CountC(ChangeRing(fp,GF(l^2)), GF(l^2));
    e2 := (a1^2 - c2) div 2;
    return true, a1, e2 - 2*l;   // <s,n> with this l's sign of trace
end function;

primes := [7,17,23,31,41,113,137,193];   // the good CRT set
Ms := AssociativeArray(); used := [];
for l in primes do
    ok, s, n := ctrlSN(l);
    if not ok then LOG(Sprintf("l=%o: control has bad reduction, skip", l)); continue; end if;
    Ml := MatchEFfast(l, s, n);
    inM := <Integers()!(GF(l)!e0), Integers()!(GF(l)!f0)> in { <x[1],x[2]> : x in Ml };
    LOG(Sprintf("l=%o <s,n>=<%o,%o> |M|=%o  true-residue-in-M? %o", l, s, n, #Ml, inM));
    Ms[l] := Ml; Append(~used, l);
end for;

// CRT reconstruct over `used`
M := &*used; ZM := Integers(M); H := Isqrt(M div 2);
LOG(Sprintf("modulus over %o primes ~%o digits, reach ~%o", #used, #IntegerToString(M), H));
np := #used; Msl := [Ms[l]:l in used]; sizes := [#m:m in Msl]; total := &*sizes;
co := [Integers()|]; for j in [1..np] do l:=used[j]; Mj:=M div l; Append(~co, Mj*(InverseMod(Mj mod l,l))); end for;
hit := false;
for cc in [0..total-1] do
    t:=cc; Es:=0; Fs:=0;
    for j in [1..np] do idx:=(t mod sizes[j])+1; t:=t div sizes[j]; pr:=Msl[j][idx];
        Es+:=pr[1]*co[j]; Fs+:=pr[2]*co[j]; end for;
    be,qe := RationalReconstruction(ZM!Es); if not be then continue; end if;
    bf,qf := RationalReconstruction(ZM!Fs); if not bf then continue; end if;
    if qe eq e0 and qf eq f0 then LOG(Sprintf("RECOVERED (e0,f0) at combo %o of %o", cc, total)); hit:=true; break; end if;
end for;
if not hit then LOG("FAILED to recover control point (check builder/reach)"); end if;
exit;
