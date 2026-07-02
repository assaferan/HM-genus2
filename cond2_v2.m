// cond2_v2.m -- conductor of Curve A at p2|2, using a MONOGENIC field (K.1 generates O_K)
// so RegularModel is available. Also try to reduce the (huge) model first.
//
// Run: magma cond2_v2.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;

Qx<X> := PolynomialRing(Rationals());
K<t> := NumberField(X^2 - X - 1);      // O_K = Z[t], t = golden ratio; s5 = 2t-1
OK := Integers(K);
s5 := 2*t - 1;
assert s5^2 eq 5;
assert Order([OK.1, OK.2]) eq OK and OK.2 eq t;  // monogenic check (t generates)

LOGF := "cond2_v2_out.txt";
PrintFile(LOGF, "# cond2_v2: conductor at 2 via monogenic field" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K);
Cmin := MinimalTwist(C, K);
Cint := IntegralModelOf(Cmin);
p2 := Factorization(2*OK)[1][1];
LOG(Sprintf("p2 norm %o; K.1 generates O_K: %o", Norm(p2), OK.2 eq t));

// try to reduce the model
Cwork := Cint;
try
    Cred := ReducedModel(Cint);
    Cwork := Cred;
    f,h := HyperellipticPolynomials(Cred);
    LOG(Sprintf("ReducedModel OK; deg f = %o, coeff heights ~ %o digits",
        Degree(f), Max([#Sprint(Numerator(c)) : c in Coefficients(f) | c ne 0])));
catch e
    LOG(Sprintf("ReducedModel unavailable/failed: %o", e`Object));
end try;

// Route A: Conductor
try
    c := Conductor(Cwork, p2);
    LOG(Sprintf("Conductor(C, p2) = %o  <-- a_2(A)", c));
catch e
    LOG(Sprintf("Conductor FAILED: %o", e`Object));
end try;

// Route B: RegularModel
try
    RM := RegularModel(Cwork, p2);
    LOG("RegularModel SUCCESS");
    try
        comps := Components(RM);
        LOG(Sprintf("  #components = %o", #comps));
    catch ee; LOG("  (could not list components)"); end try;
catch e
    LOG(Sprintf("RegularModel FAILED: %o", e`Object));
end try;

LOG("DONE.");
exit;
