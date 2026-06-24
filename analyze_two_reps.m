// analyze_two_reps.m -- Understand the two 2-dim mod-11 pieces rho_1, rho_2 of A[11].
//
// A[11] = rho_1 (+) rho_2 over K, both det = chi_11 (verified: clean-prime quadratic has
// constant N(P), linear pair has product N(P)). Then with h_P = x^4+A x^3+B x^2+C x+D the
// Frobenius quartic, the trace pair satisfies
//     {tr rho_1, tr rho_2} = roots of  y^2 + A y + (B - 2 N(P))  mod 11.
// (since t1+t2 = -A, t1 t2 = B - 2N.)  Computable at EVERY prime.
//
// Tests:
//  (1) twist hypothesis rho_2 = rho_1 (x) eps (eps quadratic): holds iff at every prime
//      A ≡ 0  OR  disc ≡ 0, where disc = A^2 - 4(B-2N).  (t2 = +-t1.)
//  (2) tabulate surface-trace (= -A) zero-rate and disc-zero-rate.
//
// Run: magma analyze_two_reps.m   (needs mod11_frobdata.txt)

SetColumns(0);
load "mod11_frobdata.txt";
F11 := GF(11);

nA0 := 0;        // primes with surface trace A ≡ 0
nDisc0 := 0;     // primes with disc ≡ 0  (t1 = t2)
nBoth := 0;      // A≡0 and disc≡0
nViol := 0;      // VIOLATIONS of twist hypothesis: A≢0 AND disc≢0
ntot := 0;
viol_examples := [];

for row in frobdata do
    Np := row[1]; p := row[2]; cs := row[4];
    A := F11 ! cs[2];           // x^3 coeff = c1
    B := F11 ! cs[3];           // x^2 coeff = c2
    N := F11 ! Np;
    disc := A^2 - 4*(B - 2*N);
    a0 := (A eq 0); d0 := (disc eq 0);
    ntot +:= 1;
    if a0 then nA0 +:= 1; end if;
    if d0 then nDisc0 +:= 1; end if;
    if a0 and d0 then nBoth +:= 1; end if;
    if not a0 and not d0 then
        nViol +:= 1;
        if #viol_examples lt 8 then Append(~viol_examples, <Np,p,Integers()!A,Integers()!disc>); end if;
    end if;
end for;

printf "Total primes: %o\n", ntot;
printf "surface trace A ≡ 0 : %o  (%.1o%%)\n", nA0, 100.0*nA0/ntot;
printf "disc ≡ 0 (t1=t2)    : %o  (%.1o%%)\n", nDisc0, 100.0*nDisc0/ntot;
printf "both                : %o\n", nBoth;
printf "\nTWIST HYPOTHESIS rho_2 = rho_1 (x) eps:  requires (A≡0 OR disc≡0) at every prime.\n";
printf "  VIOLATIONS (A≢0 and disc≢0): %o / %o\n", nViol, ntot;
if nViol eq 0 then
    printf "  => HOLDS. The two pieces are quadratic twists; one rep up to twist.\n";
else
    printf "  => FAILS. Examples <Np,p,A,disc>: %o\n", viol_examples;
    printf "  => rho_1, rho_2 are genuinely independent 2-dim reps (not a single twist class).\n";
end if;
