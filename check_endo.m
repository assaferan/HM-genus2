// check_endo.m -- Geometric endomorphism algebra of Curve A (Q(sqrt5) (1,11) curve).
// The mod-11 Frobenius factors only as [1,1,2] or [1,1,1,1] (never [2,2]/[4]/[1,3]),
// which would be explained by real multiplication with 11 split. Decide RM vs End=Z.
//
// Run: magma check_endo.m

AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
SetColumns(0);

k<s5> := QuadraticField(5);
QI := [k|
    1,
    1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895)
];
C := Genus2CurveFromIgusa(QI, k);
printf "Curve A reconstructed over Q(sqrt5).\n";

// High-level CHIMP geometric endomorphism algebra (handles embeddings, descent).
printf "Computing geometric endomorphism algebra (this can take a bit)...\n";
printf "\nGeometric End algebra: %o\n", HeuristicEndomorphismDescription(C : Geometric := true);
printf "End algebra over Q(sqrt5): %o\n", HeuristicEndomorphismDescription(C : Geometric := false);
printf "Endomorphism field of definition: %o\n", HeuristicEndomorphismFieldOfDefinition(C);
