// probe_dim.m -- find the feasibility boundary: how costly is building /
// decomposing HilbertCuspForms([2,2]) as the level norm grows?  Also test
// whether Dimension alone is cheap (mass formula) at huge norms.
SetColumns(0);
F<w> := QuadraticField(2); OF := Integers(F);
P2 := Factorization(2*OF)[1][1];   // ramified, norm 2
P3 := Factorization(3*OF)[1][1];   // 3 in Q(sqrt2): norm?
P5 := Factorization(5*OF)[1][1];   // inert, norm 25
printf "NormP2=%o NormP3=%o NormP5=%o\n", Norm(P2), Norm(P3), Norm(P5);

// mid-range level: serre 144 = 2^4*3^2 -> base * P5^2.  Build base of norm 144.
// 144 = 2^4 * 3^2.  In Q(sqrt2): 2 ramified (P2 norm2) so 2^4 -> P2^8? norm 2^8=256 no.
// Just use IdealsUpTo to grab a norm-144 ideal is overkill; construct directly:
// norm 144: need 2-part norm 16 = P2^? Norm(P2)=2 -> P2^4 (norm 16); 3-part norm 9.
f3 := Factorization(3*OF); // check 3 split/inert
printf "3 factors: %o\n", [<Norm(x[1]),x[2]> : x in f3];

for tgt in [90000, 160000, 1440000] do
    // build a level ideal of norm tgt supported on 2,3,5 (structure not critical for a timing probe)
    // tgt = base_norm * 625 where base_norm = tgt/625
    b := tgt div 625;
    // base of norm b from P2 (norm2) and P3s and ...; simplest: use P2-power * P3-power * P5^2
    // Instead: just time building at a level ideal whose norm we control loosely.
    lev := P5^2;
    n := Norm(lev);
    i2 := 0;
    while n*2 le tgt do lev := lev*P2; n := Norm(lev); i2 +:= 1; end while;
    printf "\n--- target ~%o, using level norm %o ---\n", tgt, Norm(lev);
    t0 := Cputime();
    M := HilbertCuspForms(F, lev, [2,2]);
    dm := Dimension(M);
    printf "  Dimension = %o   (%.1o s to build+dim)\n", dm, Cputime(t0);
    if Norm(lev) le 120000 then
        t1 := Cputime();
        D := NewformDecomposition(NewSubspace(M));
        printf "  NewformDecomposition: %o orbits  (%.1o s)\n", #D, Cputime(t1);
    else
        printf "  (skipping decomposition at this size)\n";
    end if;
end for;
exit;
