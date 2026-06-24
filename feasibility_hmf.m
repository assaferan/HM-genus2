// feasibility_hmf.m -- Scaling study for the Hilbert cusp forms computation over Q(sqrt5)
// at parallel weight 2, prime level. Times space setup, Dimension(NewSubspace), and one
// Hecke operator, for split primes of increasing norm, to extrapolate to level norm 66179.
//
// Run: magma feasibility_hmf.m

SetColumns(0);
K<w> := QuadraticField(5);
OK := Integers(K);

// auxiliary small prime for the Hecke-operator timing (coprime to all test levels)
qaux := Factorization(19*OK)[1][1];   // norm 19, split

// test levels: split primes (p = +-1 mod 5) of increasing norm
testnorms := [31, 89, 251, 1009, 3001, 5009];

printf "%-8o %-8o %-10o %-10o %-10o %o\n", "N(P)", "dim_new", "setup(s)", "dim(s)", "T_q(s)", "mem(MB)";
printf "%o\n", "-"^66;

for pn in testnorms do
    P := Factorization(pn*OK)[1][1];
    assert Norm(P) eq pn;
    t0 := Cputime();
    M := HilbertCuspForms(K, P, [2,2]);
    Mnew := NewSubspace(M);
    tsetup := Cputime(t0);
    t1 := Cputime();
    dn := Dimension(Mnew);
    tdim := Cputime(t1);
    t2 := Cputime();
    T := HeckeOperator(Mnew, qaux);
    tT := Cputime(t2);
    mem := GetMemoryUsage() / 1024.0 / 1024.0;
    printf "%-8o %-8o %-10o %-10o %-10o %.1o\n",
        pn, dn, RealField(4)!tsetup, RealField(4)!tdim, RealField(4)!tT, mem;
end for;

printf "\nExtrapolate the (setup+dim+Hecke) time and dimension to N(P)=66179.\n";
printf "(Full job also needs ~30-60 Hecke operators T_q reduced mod 11 to pin the eigensystem.)\n";
