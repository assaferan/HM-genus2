// goal2_mod5_qcurve2.m -- as qcurve.m but search weight-2 newforms with a
// nontrivial QUADRATIC nebentypus (the building-block descent carries the quadratic
// character of Q(sqrt2)).  Levels N=2^a*5^b, all quadratic chars mod N.
SetColumns(0);
tg := [ <7,-2,-2>,<17,6,-3>,<23,14,46>,<31,6,6>,<41,-8,-11>,<47,0,-48>,<71,4,-8>,<73,-8,13> ];
levels := [];
for a in [1..8] do for b in [0..4] do
    N := 2^a * 5^b; if N ge 8 and N le 4000 then Append(~levels, N); end if;
end for; end for;
levels := Sort(Setseq(Seqset(levels)));

LOGF := "goal2_mod5_qcurve2_out.txt";
PrintFile(LOGF, "# weight-2 newforms w/ quadratic nebentypus, Hecke field Q(sqrt3), matching B_f" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

for N in levels do
    G := FullDirichletGroup(N);
    chis := [ chi : chi in Elements(G) | Order(chi) le 2 ];
    for chi in chis do
        M := ModularForms(chi, 2);
        S := CuspidalSubspace(M);
        if Dimension(S) eq 0 then continue; end if;
        NF := Newforms(S);
        for orb in NF do
            f := orb[1]; K := BaseRing(Parent(f));
            if Type(K) eq FldRat or Degree(K) ne 2 then continue; end if;
            if Discriminant(Integers(K)) ne 12 then continue; end if;   // Q(sqrt3)
            ok := true; nt := 0;
            for t in tg do
                l := t[1]; if N mod l eq 0 then continue; end if;
                al := Coefficient(f, l);
                if (Integers()!Trace(K!al) ne t[2]) or (Integers()!Norm(K!al) ne t[3]) then ok := false; break; end if;
                nt +:= 1;
            end for;
            if ok and nt ge 4 then
                LOG(Sprintf(">>> MATCH level %o, char order %o cond %o, Hecke %o, %o primes",
                    N, Order(chi), Conductor(chi), DefiningPolynomial(K), nt));
                LOG(Sprintf("    a7=%o a17=%o a23=%o", Coefficient(f,7),Coefficient(f,17),Coefficient(f,23)));
            end if;
        end for;
    end for;
    LOG(Sprintf("  level %o done", N));
end for;
LOG("done");
exit;
