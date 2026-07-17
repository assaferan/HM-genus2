// goal2_mod5_extend.m -- extend B_f's fingerprint to more SPLIT primes of F=Q(sqrt2)
// (l = ell prime, ell split <=> ell = +-1 mod 8), recording <ell, Tr(a_P), Nm(a_P)>.
// These feed the reconstruction: we want extra LARGE primes with SMALL |M_l| to raise
// the CRT modulus (height reach ~ sqrt(prod l)) without exploding prod|M_l|.
// Also builds M_l for each and reports |M_l|, appending the small ones to mcache.m.
//
// Args: LO, HI (prime range, default 200..460).
// Run: magma LO:=200 HI:=460 goal2_mod5_extend.m
SetColumns(0);
if assigned LO then LO := StringToInteger(LO); else LO := 200; end if;
if assigned HI then HI := StringToInteger(HI); else HI := 460; end if;

F<w> := QuadraticField(2); OF := Integers(F);
P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1];
nlev := P2^3 * P5^2;
M := HilbertCuspForms(F, nlev, [2,2]);
Dc := NewformDecomposition(NewSubspace(M));
ef := Eigenform(Dc[6]);
K := HeckeEigenvalueField(Dc[6]);

LOGF := "goal2_mod5_extend_out.txt";
PrintFile(LOGF, Sprintf("# extend fingerprint + |M_l|, primes %o..%o", LO, HI) : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

// fingerprint <s,n> at the new split primes
function fingerprint(l)
    // l splits: take one prime above; a_P in K, s=Tr,n=Nm
    P := Factorization(l*OF)[1][1];
    aP := HeckeEigenvalue(ef, P);
    return Trace(K!aP), Norm(K!aP);
end function;

// --- H_12 match-set builder (same as mbuild) ---
function IGtuple(e,f)
    A1:=f-1; A:=-(f^4+15*e*f+9*e)/3; B1:=(2*f^3-2*f^2-3*f+3*e+3)/3;
    B:=(2*f^6-63*e*f^3-81*e*f^2-54*e^2)/27; B2:=e^3;
    return [ -24*B1/A1, -12*A, 96*(A/A1)*B1-36*B, -4*A1*B2 ];
end function;
badl := func< e,f | e*f*(f^2-1)*(e-f^2)*(8*f^4-9*e*f^2-3*e)*(f^6-f^4-18*e*f^2+27*e^2+16*e) >;
function CountC(fp, Fq)
    N := 0;
    for x in Fq do v := Evaluate(fp, x);
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

CACHE := "mcache_ext.m";
PrintFile(CACHE, "// extra cached M_l (small |M|, large l)" : Overwrite := true);
PrintFile(CACHE, "Mext := AssociativeArray();");

CAP := 20;   // only keep primes with |M| <= CAP (cheap combinatorics)
t0 := Cputime();
for l in PrimesInInterval(LO, HI) do
    if (l mod 8 ne 1) and (l mod 8 ne 7) then continue; end if;  // split in Q(sqrt2)
    s, n := fingerprint(l);
    tl := Cputime();
    Ml := MatchEFfast(l, s, n);
    keep := #Ml le CAP and #Ml gt 0;
    LOG(Sprintf("l=%o <s,n>=<%o,%o> |M|=%o %o [%.1o s]", l, s, n, #Ml, keep select "KEEP" else "", Cputime(tl)));
    if keep then
        entries := [ Sprintf("<%o,%o>", x[1], x[2]) : x in Ml ];
        PrintFile(CACHE, Sprintf("Mext[%o] := [ %o ];", l, Join(entries, ", ")));
    end if;
end for;
LOG(Sprintf("done in %.1o s", Cputime(t0)));
exit;
