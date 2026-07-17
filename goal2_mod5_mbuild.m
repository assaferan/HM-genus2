// goal2_mod5_mbuild.m -- build M_l (H_12 (e,f) matching B_f's L-poly up to twist) for
// ALL split primes in the fingerprint, report |M_l|, and CACHE each M_l to a loadable
// data file (mcache.m).  Reused by the solver so we never rebuild.
SetColumns(0);
function IGtuple(e,f)
    A1:=f-1; A:=-(f^4+15*e*f+9*e)/3; B1:=(2*f^3-2*f^2-3*f+3*e+3)/3;
    B:=(2*f^6-63*e*f^3-81*e*f^2-54*e^2)/27; B2:=e^3;
    return [ -24*B1/A1, -12*A, 96*(A/A1)*B1-36*B, -4*A1*B2 ];
end function;
badl := func< e,f | e*f*(f^2-1)*(e-f^2)*(8*f^4-9*e*f^2-3*e)*(f^6-f^4-18*e*f^2+27*e^2+16*e) >;

// B_f fingerprint <s,n> at split primes (L = l^2 T^4 - l s T^3 + (2l+n) T^2 - s T + 1).
target := AssociativeArray();
for t in [ <7,-2,-2>,<17,6,-3>,<23,14,46>,<31,6,6>,<41,-8,-11>,<47,0,-48>,<71,4,-8>,
           <73,-8,13>,<79,2,-2>,<89,-6,-3>,<97,-4,-8>,<103,16,16>,<113,32,253>,
           <127,18,-66>,<137,-26,121>,<151,14,-98>,<167,-12,-72>,<191,6,-138>,
           <193,34,181>,<199,8,-176> ] do target[t[1]]:=<t[2],t[3]>; end for;

function CountC(fp, Fq)
    N := 0;
    for x in Fq do
        v := Evaluate(fp, x);
        if v eq 0 then N +:= 1; elif IsSquare(v) then N +:= 2; end if;
    end for;
    d := Degree(fp);
    if IsEven(d) then if IsSquare(LeadingCoefficient(fp)) then N +:= 2; end if;
    else N +:= 1; end if;
    return N;
end function;

function MatchEFfast(l)
    Fl := GF(l); Fl2 := GF(l^2); sn := target[l]; s := sn[1]; n := sn[2];
    M := [];
    for ev in Fl do for fv in Fl do
        if badl(ev,fv) eq 0 then continue; end if;
        ic := IGtuple(ev,fv); if ic[4] eq 0 then continue; end if;
        try C := HyperellipticCurveFromIgusaClebsch(ic); catch e continue; end try;
        fp := HyperellipticPolynomials(C);
        if Degree(fp) lt 5 then continue; end if;
        a1 := l + 1 - CountC(fp, Fl);
        if AbsoluteValue(a1) ne AbsoluteValue(s) then continue; end if;
        c2 := l^2 + 1 - CountC(ChangeRing(fp, Fl2), Fl2);
        e2 := (a1^2 - c2) div 2;
        if e2 - 2*l eq n then Append(~M, <Integers()!ev, Integers()!fv>); end if;
    end for; end for;
    return M;
end function;

LOGF := "goal2_mod5_mbuild_out.txt";
PrintFile(LOGF, "# M_l build for all fingerprint primes" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

primes := [7,17,23,31,41,47,71,73,79,89,97,103,113,127,137,151,167,191,193,199];
CACHE := "mcache.m";
PrintFile(CACHE, "// cached M_l sets: Mcache[l] = list of <e,f> in Z^2 (reps mod l)" : Overwrite := true);
PrintFile(CACHE, "Mcache := AssociativeArray();");

t0 := Cputime();
for l in primes do
    tl := Cputime();
    M := MatchEFfast(l);
    LOG(Sprintf("l=%o |M|=%o  [%.1o s this, %.1o s total]", l, #M, Cputime(tl), Cputime(t0)));
    // write to cache
    entries := [ Sprintf("<%o,%o>", x[1], x[2]) : x in M ];
    PrintFile(CACHE, Sprintf("Mcache[%o] := [ %o ];", l, Join(entries, ", ")));
end for;
LOG(Sprintf("done in %.1o s", Cputime(t0)));
exit;
