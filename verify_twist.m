// verify_twist.m -- test whether the level-Pbad form f realizes rho_i up to the twist chi_g,
// where g generates the (squarefree) residual content of the model = the quadratic character
// ramified at {499,13711,88301,231481, huge primes} that the reconstruction introduced.
// If chi_g(P) matches the recorded sign eps(P), then f IS the modular form for rho_i (twisted),
// and A is a wrong-twist model: the true level is Pbad*(2)^a (no extra primes).
//
// Run: magma verify_twist.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
Qx<X> := PolynomialRing(Rationals());
K<t> := NumberField(X^2 - X - 1); OK := Integers(K); s5 := 2*t - 1;
PK<x> := PolynomialRing(K);

LOGF := "verify_twist_out.txt";
PrintFile(LOGF, "# verify_twist" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K); Cint := IntegralModelOf(Cmin);
f, h := HyperellipticPolynomials(Cint);
F0 := PK![OK!c : c in Coefficients(f + (Parent(f)!h)^2/4)];
cont := &+[ideal<OK|c> : c in Coefficients(F0) | c ne 0];
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(cont)]);
_, alpha := IsPrincipal(Msq);
F1 := PK![ (OK!c)/alpha^2 : c in Coefficients(F0) ];
cont1 := &+[ideal<OK|c> : c in Coefficients(F1) | c ne 0];
isp, g := IsPrincipal(cont1);
LOG(Sprintf("residual content principal: %o, factorization: %o", isp,
    [<Norm(pe[1]), pe[2]> : pe in Factorization(cont1)]));

function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, t - rid>;   // NB: this K uses generator t (golden), rid was red(s5); adjust below
end function;

// chi_g(P) = +1 if g is a square in the residue field at P, else -1
function ChiG(P)
    Fq, red := ResidueClassField(P);
    v := red(OK!g);
    if v eq 0 then return 0; end if;
    return IsSquare(v) select 1 else -1;
end function;

// load signs: lines "N p rid  a=.. tTr=.. tNm=..  eps=+/-/0"; rid is red(s5) in the OLD field.
// Our K here is x^2-x-1 with s5=2t-1, so red(s5)=2*rid_t-1... but the signs file rid is red(s5)
// in QuadraticField(5). To rebuild P from (p, rid_s5): P = <p, s5 - rid_s5> = <p, 2t-1-rid_s5>.
lines := Split(Read("refine_serre5_signs.txt"), "\n");
nmatch := 0; nanti := 0; ntot := 0; nzero := 0;
for ln in lines do
    if #ln eq 0 or ln[1] eq "#" then continue; end if;
    fs := [w : w in Split(ln, " ") | #w gt 0];
    p := StringToInteger(fs[2]); rid := StringToInteger(fs[3]);
    e := fs[#fs]; bit := e[#e];
    if bit ne "+" and bit ne "-" then continue; end if;
    if rid eq -1 then P := ideal<OK|p>; else P := ideal<OK | p, 2*t-1 - rid>; end if;
    if Norm(P) notin {p, p^2} then continue; end if;
    cg := ChiG(P);
    if cg eq 0 then nzero +:= 1; continue; end if;
    epsval := (bit eq "+") select 1 else -1;
    ntot +:= 1;
    if cg eq epsval then nmatch +:= 1; else nanti +:= 1; end if;
end for;
LOG(Sprintf("chi_g vs eps over %o unambiguous primes: match %o, anti-match %o (zeros %o)",
    ntot, nmatch, nanti, nzero));
if nmatch eq ntot or nanti eq ntot then
    LOG("*** chi_g REPRODUCES the sign eps EXACTLY: f realizes rho_i up to the twist chi_g. ***");
    LOG("*** => A_recon is a wrong twist; true level is Pbad*(2)^a (no 499/13711/... primes). ***");
else
    LOG("chi_g does not match eps; the twist is something else.");
end if;
LOG("DONE.");
exit;
