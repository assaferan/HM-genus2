// lseries_probe.m -- can we build the L-series of C/K and determine the conductor at 2
// from the functional equation? (K=Q(sqrt5), 2 inert; a_2(W)=2*a_{p2}(A).)
//
// Run: magma lseries_probe.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;

Qx<X> := PolynomialRing(Rationals());
K<t> := NumberField(X^2 - X - 1); OK := Integers(K); s5 := 2*t - 1;

LOGF := "lseries_probe_out.txt";
PrintFile(LOGF, "# lseries_probe" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K);
Cmin := MinimalTwist(C, K);
Cint := IntegralModelOf(Cmin);
LOG("curve built over K");

// Try to build the L-series of the curve over K
try
    L := LSeries(Cint);
    LOG("LSeries(C) built");
    try LOG(Sprintf("  Conductor(L) = %o", Conductor(L))); catch e; LOG(Sprintf("  Conductor(L) failed: %o", e`Object)); end try;
    try LOG(Sprintf("  Degree(L) = %o", Degree(L))); catch e; end try;
    try
        cfe := CheckFunctionalEquation(L);
        LOG(Sprintf("  CheckFunctionalEquation = %o", cfe));
    catch e; LOG(Sprintf("  CFE failed: %o", e`Object)); end try;
catch e
    LOG(Sprintf("LSeries(C) FAILED: %o", e`Object));
end try;

// Also try LSeries of the Jacobian
try
    J := Jacobian(Cint);
    LJ := LSeries(J);
    LOG("LSeries(Jacobian) built");
    try LOG(Sprintf("  Conductor = %o", Conductor(LJ))); catch e; end try;
catch e
    LOG(Sprintf("LSeries(Jacobian) FAILED: %o", e`Object));
end try;

LOG("DONE.");
exit;
