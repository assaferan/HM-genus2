// verify_twist.m -- confirm eps = chi_{K(sqrt d)}, d = 66179*(5-sqrt5)/2, against eps_data.
SetColumns(0);
K<w> := QuadraticField(5); OK := Integers(K);
load "eps_data.txt"; load "mod11_trnm.txt"; F11 := GF(11);
function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if; return ideal<OK | p, w - rid>; end function;
degen := {<r[2],r[3]> : r in trnm | (F11!r[4])^2 - 4*(F11!r[5]) eq 0};

d := 66179*(5-w)/2;
printf "d = %o\n", d;
printf "  v_(sqrt5)(d) = %o (odd => ramified at 5)\n", Valuation(d, Factorization(5*OK)[1][1]);
for f in Factorization(66179*OK) do
    printf "  v_(66179-prime)(d) = %o\n", Valuation(d, f[1]); end for;
printf "  v_(2)(d) = %o (0 => unramified at 2)\n\n", Valuation(d, Factorization(2*OK)[1][1]);

nND := 0; ok := 0; misfit := [];
for e in epsraw do
    if <e[1],e[2]> in degen then continue; end if;
    P := PrimeFromTag(e[1],e[2]); Fp, red := ResidueClassField(P);
    dm := red(OK!d); if dm eq 0 then continue; end if;
    sp := IsSquare(dm) select 1 else -1;
    nND +:= 1;
    if sp eq e[3] then ok +:= 1; else Append(~misfit, <e[1],e[2]>); end if;
end for;
printf "chi_{K(sqrt d)} matches %o / %o non-degenerate primes\n", ok, nND;
printf "misfit prime(s): %o\n", misfit;
