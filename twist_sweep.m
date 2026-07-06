// twist_sweep.m -- does a quadratic twist chi_d (ramified only at 2) lower a_2(rho_i) into a
// computable level? For each 2-ramified quadratic chi_d, match A_true (x) chi_d (fingerprint =
// (chi_g*chi_d)-twist of mod11_trnm) via untwisted V+ at levels Pbad*(2)^{0,1,2}. Hecke ops
// computed once per level; the chi_d sweep is cheap (sign flips + incremental kernel intersect).
//
// Run: magma twist_sweep.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf, BadPrimesFromInvariants;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K); F11 := GF(11);
LOGF := "twist_sweep_out.txt"; PrintFile(LOGF, "# twist_sweep" : Overwrite := true);
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

// enumerate 2-ramified quadratic characters via the ray class group mod p2^6 * oo
R, mR := RayClassGroup(p2^6, [1,2]);
Q, qR := quo<R | [2*R.i : i in [1..Ngens(R)]]>;
rk := Ngens(Q);
LOG(Sprintf("rank R/2R at p2^6*oo = %o => %o quadratic characters (incl trivial)", rk, 2^rk));
// each functional f in F_2^rk gives chi_d(P) = (-1)^<f, coords(P)>
funcs := [ [GF(2)!((i div 2^(j-1)) mod 2) : j in [1..rk]] : i in [0..2^rk-1] ];

load "mod11_trnm.txt";
targets := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    P := PrimeFromTag(row[2], row[3]); cg := ChiG(P);
    if cg eq 0 then continue; end if;
    coords := [GF(2)!c : c in Eltseq(qR(P @@ mR))];
    Append(~targets, <P, cg, F11!row[4], F11!row[5], row[1], coords>);   // <P, chi_g, tTr, tNm, N, chi_d-coords>
end for;
Sort(~targets, func<a,b | a[5]-b[5]>);

procedure SweepLevel(Level, NPRIMES)
    LOG(Sprintf("==== level norm %o ====", Norm(Level)));
    M := HilbertCuspForms(K, Level, [2,2]); d := Dimension(M); Id := IdentityMatrix(F11,d);
    LOG(Sprintf("  dim %o; computing Hecke ops...", d));
    Ts := []; tgs := [];
    for tg in targets do
        if Norm(tg[1]+Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        Append(~Ts, T); Append(~tgs, tg);
        if #Ts ge NPRIMES then break; end if;
    end for;
    LOG(Sprintf("  %o Hecke ops ready; sweeping %o twists...", #Ts, #funcs));
    for fi in [1..#funcs] do
        fnc := funcs[fi];
        Vp := VectorSpace(F11, d); alive := true;
        for j in [1..#Ts] do
            tg := tgs[j]; T := Ts[j];
            chid := IsEven(Integers()!(&+[Integers()| fnc[k]*Integers()!tg[6][k] : k in [1..rk]])) select F11!1 else F11!(-1);
            sgn := tg[2]*chid;                     // chi_g(P)*chi_d(P)
            Vp := Vp meet Kernel(T*T - sgn*tg[3]*T + tg[4]*Id);
            if Dimension(Vp) eq 0 then alive := false; break; end if;
        end for;
        if alive and Dimension(Vp) gt 0 then
            LOG(Sprintf("  *** twist #%o: V+ = %o (SURVIVES) -- a_2 drops here! ***", fi, Dimension(Vp)));
        end if;
    end for;
    LOG("  (no surviving twist unless noted above)");
end procedure;

SweepLevel(Pbad, 12);
SweepLevel(Pbad*p2, 12);
SweepLevel(Pbad*p2^2, 10);
LOG("DONE."); exit;
