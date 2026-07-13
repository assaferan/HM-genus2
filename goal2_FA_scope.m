// goal2_FA_scope.m -- scope the F_A resolvent computation for idx 33.
// (1) Characterize the target: image of sigma_A and the degree-50 PSL_2(F_49)
//     projective resolvent.  (2) Probe whether the genus-2 7-torsion is computable.
SetColumns(0);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
C := HyperellipticCurve(-10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9);
LOGF := "goal2_FA_scope_out.txt";
PrintFile(LOGF, "# F_A scoping (idx 33)" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

// ---- (1) target image group and projective resolvent (orders computed analytically) ----
q := 49;
ordGL  := (q^2-1)*(q^2-q);            // |GL_2(F_49)|
ordSL  := q*(q^2-1);                  // |SL_2(F_49)|
imord  := ordGL * 6 div 48;           // {det in F_7^*}: preimage of order-6 in F_49^* (order 48)
ordPSL := ordSL div 2;                // |PSL_2(F_49)|
LOG("=== target: image and projective resolvent ===");
LOG(Sprintf("  |GL_2(F_49)| = %o,  |SL_2(F_49)| = %o", ordGL, ordSL));
LOG(Sprintf("  im(sigma_A) = { g in GL_2(F_49) : det g in F_7^* },  order %o", imord));
LOG(Sprintf("  projective image = PSL_2(F_49),  order %o", ordPSL));
LOG(Sprintf("  PSL_2(F_49) acts on P^1(F_49) = 50 points (transitive), stabilizer order %o", ordPSL div 50));
LOG("  => F_A projective resolvent: degree 50 over F=Q(sqrt10), Galois group PSL_2(F_49),");
LOG("     ramified only at {2,3,5,7}; over Q the tau-swap extends it by Frob_{F49/F7}.");
LOG(Sprintf("  (full 7-torsion field has degree |image A[7]| = %o over Q.)", 2*imord));

// ---- (2) probe genus-2 7-torsion tooling ----
LOG("\n=== probe: 7-torsion tooling ===");
J := Jacobian(C);
LOG("  Jacobian built.");
t0 := Cputime();
try
    K := KummerSurface(J);
    LOG(Sprintf("  KummerSurface built (%.1o s).", Cputime(t0)));
catch e; LOG(Sprintf("  KummerSurface failed: %o", e`Object)); end try;
try
    dp := DivisionPolynomial(J, 7);
    LOG(Sprintf("  DivisionPolynomial(J,7): degree %o", Degree(dp)));
catch e; LOG(Sprintf("  DivisionPolynomial(J,7) not available: %o", e`Object)); end try;
LOG("DONE.");
exit;
