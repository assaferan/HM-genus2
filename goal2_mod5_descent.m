// goal2_mod5_descent.m -- Quer descent: find the classical weight-2 newform g/Q
// (ANY Hecke degree) whose TRACE FORM matches Res(B_f): Tr(a_l(g)) = s_l at primes
// l split in Q(sqrt2), = 0 at inert l.  Nebentypus is the quadratic char of Q(sqrt2)
// (chi_8).  This is the descended object; Res(B_f) ~ A_g.
SetColumns(0);
// trace targets: Tr(a_l) at good primes.  split -> s_l ; inert -> 0.
// split l of Q(sqrt2): 7,17,23,31,41,47,71,73 (s from fingerprint); inert: 3,11,13,19,29,37,43
sp := [ <7,-2>,<17,6>,<23,14>,<31,6>,<41,-8>,<47,0> ];
inert := [3,11,13,19,29,37,43];
tr := AssociativeArray();
for t in sp do tr[t[1]] := t[2]; end for;
for l in inert do tr[l] := 0; end for;
testp := Sort([ k : k in Keys(tr) ]);

levels := [];
for a in [1..8] do for b in [0..4] do
    N := 2^a*5^b; if N ge 8 and N le 3200 then Append(~levels, N); end if;
end for; end for;
levels := Sort(Setseq(Seqset(levels)));

LOGF := "goal2_mod5_descent_out.txt";
PrintFile(LOGF, "# Quer descent: classical newform g with Res(B_f)~A_g, trace form supported on split primes" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

for N in levels do
    G := FullDirichletGroup(N);
    for chi in [ c : c in Elements(G) | Order(c) le 2 ] do
        M := ModularForms(chi, 2); S := CuspidalSubspace(M);
        if Dimension(S) eq 0 then continue; end if;
        NF := Newforms(S);
        for orb in NF do
            f := orb[1];
            ok := true; nt := 0;
            for l in testp do
                if N mod l eq 0 then continue; end if;
                al := Coefficient(f, l);
                K := Parent(al);
                trv := (Type(K) eq FldRat) select (Integers()!al) else (Integers()!Trace(al));
                if trv ne tr[l] then ok := false; break; end if;
                nt +:= 1;
            end for;
            if ok and nt ge 5 then
                K := BaseRing(Parent(f));
                deg := (Type(K) eq FldRat) select 1 else Degree(K);
                LOG(Sprintf(">>> MATCH level %o, char cond %o order %o, Hecke deg %o, %o primes",
                    N, Conductor(chi), Order(chi), deg, nt));
                LOG(Sprintf("    a3=%o a7=%o a11=%o a17=%o a23=%o",
                    Coefficient(f,3),Coefficient(f,7),Coefficient(f,11),Coefficient(f,17),Coefficient(f,23)));
            end if;
        end for;
    end for;
    LOG(Sprintf("  level %o done", N));
end for;
LOG("done");
exit;
