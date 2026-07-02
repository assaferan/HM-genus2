// minimize1.m -- reduce Curve A's model over K: remove square content (isomorphism y->a*y),
// report remaining (squarefree) content and coefficient size, then try coefficient reduction.
//
// Run: magma minimize1.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;

Qx<X> := PolynomialRing(Rationals());
K<t> := NumberField(X^2 - X - 1); OK := Integers(K); s5 := 2*t - 1;

LOGF := "minimize1_out.txt";
PrintFile(LOGF, "# minimize1" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K); Cint := IntegralModelOf(Cmin);
f, h := HyperellipticPolynomials(Cint);
PK := Parent(f);
F0 := f + (PK!h)^2/4;   // y^2 = F0 (h=0 here anyway)
F0 := PK![OK!c : c in Coefficients(F0)];

function ModelSize(poly)
    return Max([#Sprint(Numerator(Norm(c))) : c in Coefficients(poly) | c ne 0]);
end function;
LOG(Sprintf("orig: deg %o, max |coeff| ~ %o digits", Degree(F0), ModelSize(F0)));

cont := &+[ideal<OK|c> : c in Coefficients(F0) | c ne 0];
fac := Factorization(cont);
LOG(Sprintf("content factorization (norm:exp): %o", [<Norm(pe[1]), pe[2]> : pe in fac]));

// square part: M = prod P^floor(e/2); alpha^2 = M^2 (principal, h(K)=1)
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in fac]);
issq, alpha := IsPrincipal(Msq);
assert issq;
LOG(Sprintf("removing alpha^2, N(alpha) = %o", Norm(alpha)));

F1 := PK![ (OK!c) / alpha^2 : c in Coefficients(F0) ];
// verify integral
assert &and[ c in OK : c in Coefficients(F1) ];
LOG(Sprintf("after square-content removal: max |coeff| ~ %o digits", ModelSize(F1)));
cont1 := &+[ideal<OK|c> : c in Coefficients(F1) | c ne 0];
LOG(Sprintf("remaining content (norm:exp): %o", [<Norm(pe[1]), pe[2]> : pe in Factorization(cont1)]));

// sanity: same curve (Igusa invariants match up to scaling)
C1 := HyperellipticCurve(F1);
IC0 := IgusaClebschInvariants(Cint);
IC1 := IgusaClebschInvariants(C1);
// compare as weighted-projective point (ratios)
function WPSeq(IC) return [IC[i]/IC[1]^(w[i] div 2) : i in [2..#IC]] where w := [2,4,6,10]; end function;
same := false;
try same := WPSeq(IC0) eq WPSeq(IC1); catch e; end try;
LOG(Sprintf("Igusa-Clebsch invariants preserved: %o", same));
LOG("F1 (reduced-content model):");
LOG(Sprintf("  %o", F1));
LOG("DONE.");
exit;
