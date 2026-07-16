// goal2_mod5_loo2.m -- FAST leave-one-out rational reconstruction of B_f's H_12
// moduli point.  Speedups vs the naive builder:
//  (1) trace-prune: compute #C(F_l) (O(l)) first; only compute #C(F_{l^2}) (O(l^2),
//      the norm) for the few (e,f) whose |trace| already matches -> avoids O(l^4).
//  (2) only the 7 small CRT primes (<=79); verify survivors at the DROPPED prime
//      (membership) + a couple fresh primes built on-the-fly (few survivors).
SetColumns(0);
function IGtuple(e,f)
    A1:=f-1; A:=-(f^4+15*e*f+9*e)/3; B1:=(2*f^3-2*f^2-3*f+3*e+3)/3;
    B:=(2*f^6-63*e*f^3-81*e*f^2-54*e^2)/27; B2:=e^3;
    return [ -24*B1/A1, -12*A, 96*(A/A1)*B1-36*B, -4*A1*B2 ];
end function;
badl := func< e,f | e*f*(f^2-1)*(e-f^2)*(8*f^4-9*e*f^2-3*e)*(f^6-f^4-18*e*f^2+27*e^2+16*e) >;
target := AssociativeArray();
for t in [ <7,-2,-2>,<17,6,-3>,<23,14,46>,<31,6,6>,<41,-8,-11>,<47,0,-48>,<71,4,-8>,
           <73,-8,13>,<79,2,-2>,<89,-6,-3>,<97,-4,-8>,<103,16,16>,<113,32,253>,
           <127,18,-66>,<137,-26,121>,<151,14,-98>,<167,-12,-72>,<191,6,-138>,
           <193,34,181>,<199,8,-176> ] do target[t[1]]:=<t[2],t[3]>; end for;

// count #C(F_q) for y^2=f(x) via Legendre (f over F_q, q a prime power)
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

// M_l with trace-pruning.  sn=(s,n): match |a1|=|s| and n exactly (up to twist).
function MatchEFfast(l)
    Fl := GF(l); Fl2 := GF(l^2); sn := target[l]; s := sn[1]; n := sn[2];
    M := [];
    for ev in Fl do for fv in Fl do
        if badl(ev,fv) eq 0 then continue; end if;
        ic := IGtuple(ev,fv); if ic[4] eq 0 then continue; end if;
        try C := HyperellipticCurveFromIgusaClebsch(ic); catch e continue; end try;
        fp := HyperellipticPolynomials(C);
        if Degree(fp) lt 5 then continue; end if;
        a1 := l + 1 - CountC(fp, Fl);            // trace
        if AbsoluteValue(a1) ne AbsoluteValue(s) then continue; end if;   // prune
        c2 := l^2 + 1 - CountC(ChangeRing(fp, Fl2), Fl2);  // sum alpha_i^2
        e2 := (a1^2 - c2) div 2;                  // 2nd elem symm
        if e2 - 2*l eq n then Append(~M, <ev,fv>); end if;   // norm match
    end for; end for;
    return M;
end function;

LOGF:="goal2_mod5_loo2_out.txt";
PrintFile(LOGF,"# fast leave-one-out reconstruction":Overwrite:=true);
procedure LOG(s) printf "%o\n",s; PrintFile(LOGF,s); end procedure;

allp := [7,23,17,31,41,73,79];
fresh := [47,89,97,103,113,127,137];  // scoring primes (all with manageable |M|)
Mset := AssociativeArray(); Msetset := AssociativeArray();
t0 := Cputime();
for l in allp cat fresh do
    Mset[l] := MatchEFfast(l);
    Msetset[l] := { x : x in Mset[l] };
    LOG(Sprintf("l=%o |M|=%o  [%.1o s]", l, #Mset[l], Cputime(t0)));
end for;
function memOK(e0,f0,l)
    if (Denominator(e0) mod l ne 0) and (Denominator(f0) mod l ne 0)
       and <GF(l)!e0,GF(l)!f0> in Msetset[l] then return true; end if;
    return false;
end function;

// fresh verification: build curve mod l for a candidate (few survivors), check L-poly
RTq<Tq> := PolynomialRing(Rationals());
Ltgt := func< l,s,n | l^2*Tq^4 - l*s*Tq^3 + (2*l+n)*Tq^2 - s*Tq + 1 >;
function checkL(e0,f0,l)
    Fl := GF(l); PF<TT> := PolynomialRing(Fl); sn := target[l];
    Lt := PF!Ltgt(l,sn[1],sn[2]); Lw := PF!Ltgt(l,-sn[1],sn[2]);
    if (Denominator(e0) mod l eq 0) or (Denominator(f0) mod l eq 0) then return false; end if;
    if badl(Fl!e0,Fl!f0) eq 0 then return false; end if;
    ic := IGtuple(Fl!e0,Fl!f0); if ic[4] eq 0 then return false; end if;
    try C := HyperellipticCurveFromIgusaClebsch([Fl|x:x in ic]); catch e return false; end try;
    return LPolynomial(C) in {Lt,Lw};
end function;

function reconstruct(used)
    Ms := [* Mset[l] : l in used *]; np:=#used;
    sizes := [ #Ms[j] : j in [1..np] ]; total := &*sizes;
    Mprod := &*used; ZM := Integers(Mprod); sols := {};
    for cc in [0..total-1] do
        t:=cc; sel:=[];
        for j in [1..np] do Append(~sel,(t mod sizes[j])+1); t:=t div sizes[j]; end for;
        be,qe := RationalReconstruction(ZM!CRT([Integers()!Ms[j][sel[j]][1]:j in [1..np]],used));
        if not be then continue; end if;
        bf,qf := RationalReconstruction(ZM!CRT([Integers()!Ms[j][sel[j]][2]:j in [1..np]],used));
        if not bf then continue; end if;
        Include(~sols,<qe,qf>);
    end for;
    return sols, Isqrt(Mprod div 2);
end function;

// Robust: reconstruct over each drop-subset, then SCORE every candidate by how many
// of the 12 fresh primes it matches (membership).  The true (e,f) scores near-max
// (missing only its badlocus primes); spurious candidates score low.  Track global best.
function trydropf(dropset)
    used := [ p : p in allp | p notin dropset ];
    if &*used lt 10^8 then return 0, {}; end if;
    sols := reconstruct(used);
    loc := 0; locbest := {};
    for ef in sols do
        sc := #[ l : l in fresh | memOK(ef[1],ef[2],l) ];
        if sc gt loc then loc := sc; locbest := {ef}; elif sc eq loc then Include(~locbest, ef); end if;
    end for;
    LOG(Sprintf("drop %o | recon %o | best fresh-score %o / %o", dropset, #sols, loc, #fresh));
    return loc, locbest;
end function;

dropsets := [ [Integers()|] ];
for i in [1..#allp] do Append(~dropsets, [allp[i]]); end for;
for i in [1..#allp] do for j in [i+1..#allp] do Append(~dropsets, [allp[i],allp[j]]); end for; end for;

bestscore := 0; bestlist := {};
for dropset in dropsets do
    loc, lb := trydropf(dropset);
    if loc gt bestscore then bestscore := loc; bestlist := lb;
    elif (loc eq bestscore) and (loc gt 0) then bestlist := bestlist join lb; end if;
end for;

LOG(Sprintf("\n=== GLOBAL BEST fresh-score %o / %o, %o candidate(s) ===", bestscore, #fresh, #bestlist));
for ef in bestlist do
    passed := [ l : l in fresh | memOK(ef[1],ef[2],l) ];
    LOG(Sprintf("  (e,f) = (%o, %o)  [passes %o]", ef[1], ef[2], passed));
end for;
exit;
