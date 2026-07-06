// verify_twist7.m -- verify the surviving quadratic twist (#7) at level Pbad over MANY primes,
// identify it, and test whether the level-Pbad form genuinely realizes rho_i up to that twist
// (consistent piece) vs the a_P^2 coincidence (jumping piece).
//
// Run: magma verify_twist7.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf, BadPrimesFromInvariants;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K); F11 := GF(11);
LOGF := "verify_twist7_out.txt"; PrintFile(LOGF, "# verify_twist7" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787), 1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089), 1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K); Cint := IntegralModelOf(Cmin);
ff,hh := HyperellipticPolynomials(Cint);
F0 := PK![OK!c : c in Coefficients(ff + (Parent(ff)!hh)^2/4)];
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(&+[ideal<OK|c> : c in Coefficients(F0)|c ne 0])]);
_, alpha := IsPrincipal(Msq);
F1 := PK![ (OK!c)/alpha^2 : c in Coefficients(F0) ];
_, g := IsPrincipal(&+[ideal<OK|c> : c in Coefficients(F1)|c ne 0]);
Pbad := [P : P in BadPrimesFromInvariants(QI, K) | Norm(P) eq 66179][1];
p2 := Factorization(2*OK)[1][1];
function PrimeFromTag(p, rid) if rid eq -1 then return ideal<OK|p>; end if; return ideal<OK|p, s5-rid>; end function;
function ChiG(P) Fp,red := ResidueClassField(P); v := red(OK!g); if v eq 0 then return 0; end if; return IsSquare(v) select F11!1 else F11!(-1); end function;

R, mR := RayClassGroup(p2^6, [1,2]);
Q, qR := quo<R | [2*R.i : i in [1..Ngens(R)]]>; rk := Ngens(Q);
f7 := [GF(2)!((6 div 2^(j-1)) mod 2) : j in [1..rk]];   // twist #7 = index 6 (0-based)
LOG(Sprintf("twist #7 functional on R/2R (rk %o): %o", rk, f7));
function ChiD7(P) coords := [GF(2)!c : c in Eltseq(qR(P @@ mR))];
  return IsEven(Integers()!(&+[Integers()| f7[k]*Integers()!coords[k] : k in [1..rk]])) select F11!1 else F11!(-1); end function;

load "mod11_trnm.txt";
targets := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    P := PrimeFromTag(row[2], row[3]); cg := ChiG(P); if cg eq 0 then continue; end if;
    Append(~targets, <P, cg*ChiD7(P), F11!row[4], F11!row[5], row[1]>);   // combined sign chi_g*chi_d7
end for;
Sort(~targets, func<a,b|a[5]-b[5]>);

// V+ over MANY primes at level Pbad
M := HilbertCuspForms(K, Pbad, [2,2]); d := Dimension(M); Id := IdentityMatrix(F11,d);
LOG(Sprintf("level Pbad dim %o; V+ with twist #7 over primes:", d));
Vp := VectorSpace(F11,d); BM := 0; used := 0;
for tg in targets do
    if Norm(tg[1]+Pbad) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
    Vp := Vp meet Kernel(T*T - tg[2]*tg[3]*T + tg[4]*Id);
    used +:= 1;
    if used mod 20 eq 0 or Dimension(Vp) eq 0 then LOG(Sprintf("  after %o primes: dim V+ = %o", used, Dimension(Vp))); end if;
    if Dimension(Vp) eq 0 then break; end if;
    if used ge 80 then break; end if;
end for;
LOG(Sprintf("=> V+ with twist #7 = %o after %o primes", Dimension(Vp), used));
if Dimension(Vp) gt 0 then
    LOG("*** GENUINE MATCH: rho_i modular by an explicit level-Pbad form, up to the quadratic twist chi_g*chi_d7. ***");
    LOG("    (the earlier a_P^2 'red herring' was mis-diagnosed: the twist is ramified at 499,13711,88301,231481,2 -- outside the conductor set tested by refine_charfit.)");
end if;

// identify chi_d7 as a quadratic field: find small d with (d) supported at p2 and chi_d = chi_d7
LOG("identifying chi_d7 (quadratic character ramified at 2):");
_, pigen := IsPrincipal(p2);
for d in [K| -1, pigen, -pigen, s5, -s5, pigen*s5, 1+s5, 1-s5, (1+s5)/2, 3+s5 ] do
    if d eq 0 then continue; end if;
    ok := true; nchk := 0;
    for tg in targets[1..30] do
        P := tg[1]; Fp,red := ResidueClassField(P); v := red(OK!d);
        if v eq 0 then continue; end if;
        cd := IsSquare(v) select F11!1 else F11!(-1);
        if cd ne ChiD7(P) then ok := false; break; end if; nchk +:= 1;
    end for;
    if ok and nchk ge 20 then LOG(Sprintf("  chi_d7 = quadratic character of K(sqrt(%o))", d)); end if;
end for;
LOG("DONE."); exit;
