// cond2_probe.m -- attempt to compute the conductor of Curve A at the inert prime p2 | 2.
// Tries several Magma routes (they are known to be fragile at p=2 over number fields).
//
// Run: magma cond2_probe.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
K<s5> := QuadraticField(5); OK := Integers(K);

LOGF := "cond2_probe_out.txt";
PrintFile(LOGF, "# cond2_probe" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K);
Cmin := MinimalTwist(C, K);
Cint := IntegralModelOf(Cmin);
p2 := Factorization(2*OK)[1][1];
LOG(Sprintf("p2 = inert prime above 2, norm %o", Norm(p2)));
f,h := HyperellipticPolynomials(Cint);
LOG(Sprintf("model: y^2 + (%o) y = %o", h, f));

// Route 1: Conductor(C, p2)
try
    c := Conductor(Cint, p2);
    LOG(Sprintf("Route1 Conductor(C,p2) = %o", c));
catch e
    LOG(Sprintf("Route1 FAILED: %o", e`Object));
end try;

// Route 2: RegularModel at p2, then read special fibre
try
    RM := RegularModel(Cint, p2);
    LOG("Route2 RegularModel: SUCCESS");
    try LOG(Sprintf("  components / special fibre info: %o", RM)); catch ee; end try;
catch e
    LOG(Sprintf("Route2 RegularModel FAILED: %o", e`Object));
end try;

// Route 3: base change to completion K_p2 and try there
try
    Kp2, mp := Completion(K, p2 : Precision := 200);
    Rp := PolynomialRing(Kp2);
    fp := Rp![mp(c) : c in Coefficients(f)];
    hp := Rp![mp(c) : c in Coefficients(h)];
    Cp := HyperellipticCurve(fp, hp);
    LOG("Route3 base change to K_p2: curve built");
    try
        RMp := RegularModel(Cp, UniformizingElement(Kp2));
        LOG("  RegularModel over completion: SUCCESS");
    catch ee
        LOG(Sprintf("  RegularModel over completion FAILED: %o", ee`Object));
    end try;
catch e
    LOG(Sprintf("Route3 FAILED: %o", e`Object));
end try;

LOG("DONE.");
exit;
