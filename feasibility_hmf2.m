// feasibility_hmf2.m -- confirm feasibility directly at/near the target level norm 66179.
SetColumns(0);
K<w> := QuadraticField(5); OK := Integers(K);

// a handful of small auxiliary primes for Hecke timing (simulate the full fingerprint job)
auxnorms := [19, 29, 31, 41, 59, 61, 71, 79];   // split primes, small norm
auxP := [ Factorization(n*OK)[1][1] : n in auxnorms ];

for pn in [20011, 66179] do
    P := Factorization(pn*OK)[1][1];
    assert Norm(P) eq pn;
    printf "=== level norm %o ===\n", pn;
    t0 := Cputime();
    M := HilbertCuspForms(K, P, [2,2]);
    Mnew := NewSubspace(M);
    dn := Dimension(Mnew);
    printf "  setup+dim: %os, dim_new = %o\n", RealField(4)!Cputime(t0), dn;
    // time a batch of Hecke operators (the full job needs ~40 of these)
    t1 := Cputime();
    nop := 0;
    for Q in auxP do
        if Norm(Q + P) ne 1 then continue; end if;  // skip if NOT coprime to level
        _ := HeckeOperator(Mnew, Q);
        nop +:= 1;
    end for;
    tops := Cputime(t1);
    printf "  %o Hecke operators: %os total, %os each\n", nop, RealField(4)!tops, RealField(4)!(tops/nop);
    printf "  => est. for 40 Hecke ops: %o s; mem %.1o MB\n",
        RealField(4)!(40*tops/nop), GetMemoryUsage()/1024.0/1024.0;
end for;
