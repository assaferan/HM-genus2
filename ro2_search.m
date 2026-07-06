// ro2_search.m -- find rho_2's form. Separate the pieces using g (rho_1 = g (x) chi):
// t_1(P) = chi(P)*a_P(g), t_2(P) = Tr_P - t_1(P). Then search rho_2 (x) (chi_g*chi_d) at level
// Pbad*(2)^2 via ker(T_P - chi_g(P)chi_d(P)*t_2(P)) over the 8 two-ramified twists chi_d.
//
// Run: magma ro2_search.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf, BadPrimesFromInvariants;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K); F11 := GF(11);
LOGF := "ro2_search_out.txt"; PrintFile(LOGF, "# ro2_search" : Overwrite := true);
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
function ChiG(P) Fp,red := ResidueClassField(P); v := red(OK!g); return (v ne 0 and IsSquare(v)) select F11!1 else F11!(-1); end function;
function ChiM1(P) Fp,red := ResidueClassField(P); return IsSquare(red(OK!(-1))) select F11!1 else F11!(-1); end function;

R, mR := RayClassGroup(p2^6, [1,2]); Q, qR := quo<R | [2*R.i : i in [1..Ngens(R)]]>; rk := Ngens(Q);
function ChiD(P, fnc) coords := [GF(2)!c : c in Eltseq(qR(P @@ mR))];
  return IsEven(Integers()!(&+[Integers()| fnc[k]*Integers()!coords[k] : k in [1..rk]])) select F11!1 else F11!(-1); end function;
funcs := [ [GF(2)!((i div 2^(j-1)) mod 2) : j in [1..rk]] : i in [0..2^rk-1] ];
f7 := [GF(2)!((6 div 2^(j-1)) mod 2) : j in [1..rk]];   // chi_{-1}

load "mod11_trnm.txt";
prs := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    P := PrimeFromTag(row[2], row[3]);
    Append(~prs, <P, F11!row[4], F11!row[5], row[1]>);   // <P, Tr, Nm, N>
end for;
Sort(~prs, func<a,b|a[4]-b[4]>);

// ---- Part 1: t_2(P) from g (level Pbad, V+ with twist #7) ----
M0 := HilbertCuspForms(K, Pbad, [2,2]); d0 := Dimension(M0); Id0 := IdentityMatrix(F11,d0);
Vp := VectorSpace(F11,d0); T0 := AssociativeArray();
for pr in prs do
    if Norm(pr[1]+Pbad) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M0, pr[1])), F11); T0[pr[4]] := T;
    chi := ChiG(pr[1])*ChiM1(pr[1]);
    Vp := Vp meet Kernel(T*T - chi*pr[2]*T + pr[3]*Id0);
    if Dimension(Vp) le 1 and #Keys(T0) ge 20 then break; end if;
end for;
LOG(Sprintf("level Pbad: V+ (twist #7) dim = %o", Dimension(Vp)));
v := Basis(Vp)[1];
t2 := AssociativeArray();
for pr in prs do
    if not IsDefined(T0, pr[4]) then continue; end if;
    w := v*T0[pr[4]]; k := 1; while v[k] eq 0 do k+:=1; end while; aPg := w[k]/v[k];
    chi := ChiG(pr[1])*ChiM1(pr[1]);
    t1 := chi*aPg; t2[pr[4]] := pr[2] - t1;   // t_2 = Tr - t_1
end for;
LOG(Sprintf("computed t_2 at %o primes", #Keys(t2)));

// ---- Part 2: search rho_2 (x) (chi_g*chi_d) at Pbad*(2)^2 ----
Level := Pbad*p2^2;
M2 := HilbertCuspForms(K, Level, [2,2]); d2 := Dimension(M2); Id2 := IdentityMatrix(F11,d2);
LOG(Sprintf("level Pbad*(2)^2 dim %o; computing Hecke ops...", d2));
Ts := []; tps := [];
for pr in prs do
    if Norm(pr[1]+Level) ne 1 or not IsDefined(t2, pr[4]) then continue; end if;
    Append(~Ts, ChangeRing(Matrix(HeckeOperator(M2, pr[1])), F11)); Append(~tps, pr);
    if #Ts ge 10 then break; end if;
end for;
LOG(Sprintf("  %o Hecke ops; sweeping %o twists for rho_2...", #Ts, #funcs));
for fi in [1..#funcs] do
    fnc := funcs[fi]; Vp2 := VectorSpace(F11,d2); alive := true;
    for j in [1..#Ts] do
        pr := tps[j]; sgn := ChiG(pr[1])*ChiD(pr[1], fnc);   // chi_g*chi_d
        ev := sgn*t2[pr[4]];                                  // target eigenvalue for rho_2 (x) twist
        Vp2 := Vp2 meet Kernel(Ts[j] - ev*Id2);
        if Dimension(Vp2) eq 0 then alive := false; break; end if;
    end for;
    if alive and Dimension(Vp2) gt 0 then
        LOG(Sprintf("  *** twist #%o: V+ = %o -- rho_2's form appears at Pbad*(2)^2! ***", fi, Dimension(Vp2)));
    end if;
end for;
LOG("(no surviving twist for rho_2 at this level unless noted)");
LOG("DONE."); exit;
