// Example / smoke test for the number-field case: generate a genus-2 curve over
// Q(sqrt2) from x = [1, sqrt2, 3, 4], using the unified Genus2Curve(x : K := ...).
//
// Run:  magma Genus2Curve.m test_quadratic.m
// (Takes a few minutes: the basin search is stochastic, ~10% of restarts succeed,
//  and the Mestre descent over a number field is not cheap.)

k<s2> := QuadraticField(2);
x := [k| 1, s2, 3, 4];

C, QI := Genus2Curve(x : K := k);     // Embedding := 1, i.e. sqrt2 -> +1.41421...

printf "Absolute Igusa invariants over Q(sqrt2):\n";
for i in [1..#QI] do printf "  I%o = %o\n", i, QI[i]; end for;
printf "Reconstructed curve over %o:\n%o\n", BaseRing(C), C;

assert BaseRing(C) cmpeq k;           // the curve descends to Q(sqrt2)
print "PASS: curve reconstructed over Q(sqrt2).";
