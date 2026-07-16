// goal2_mod5_qcurve.m -- find the classical weight-2 newform g/Q with Hecke field
// Q(sqrt3) that is the descent of the building block B_f.  At split primes l of
// Q(sqrt2), a_l(g) = a_P(f), i.e. Tr(a_l(g))=s_l, Nm(a_l(g))=n_l.  Search levels
// N = 2^a * 5^b (the bad primes) for a weight-2 newform with these eigenvalues.
SetColumns(0);
// (l, s, n) at split primes
tg := [ <7,-2,-2>,<17,6,-3>,<23,14,46>,<31,6,6>,<41,-8,-11>,<47,0,-48>,<71,4,-8>,<73,-8,13> ];

levels := [];
for a in [0..9] do for b in [0..5] do
    N := 2^a * 5^b;
    if N ge 20 and N le 40000 then Append(~levels, N); end if;
end for; end for;
levels := Sort(Setseq(Seqset(levels)));

LOGF := "goal2_mod5_qcurve_out.txt";
PrintFile(LOGF, "# search weight-2 newforms over Q with Hecke field Q(sqrt3), a_l matching B_f" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

for N in levels do
    S := CuspForms(N);
    NF := Newforms(S);
    for orb in NF do
        f := orb[1];
        K := BaseRing(Parent(f));
        if Type(K) eq FldRat then continue; end if;      // rational -> not Q(sqrt3)
        if Degree(K) ne 2 then continue; end if;
        if Discriminant(Integers(K)) ne 12 then continue; end if;  // Q(sqrt3), disc 12
        // check a_l trace/norm at split primes coprime to N
        ok := true; ntest := 0;
        for t in tg do
            l := t[1];
            if N mod l eq 0 then continue; end if;
            al := Coefficient(f, l);
            if (Integers()!Trace(K!al) ne t[2]) or (Integers()!Norm(K!al) ne t[3]) then ok := false; break; end if;
            ntest +:= 1;
        end for;
        if ok and ntest ge 4 then
            LOG(Sprintf(">>> MATCH: level %o, Hecke field %o, tested %o split primes",
                N, DefiningPolynomial(K), ntest));
            LOG(Sprintf("    a_7=%o a_17=%o a_23=%o", Coefficient(f,7), Coefficient(f,17), Coefficient(f,23)));
        end if;
    end for;
    LOG(Sprintf("  ... level %o done", N));
end for;
LOG("search complete");
exit;
