// ro2_search3.m -- search rho_2's form at level Pbad*(2)^3 (dim ~88000).
// t_2 from g (rho_1 = g (x) chi); match rho_2 (x) (chi_g*chi_d) via ker(T_P - chi_g chi_d t_2).
// Memory-efficient: store ~5 Hecke ops; restrict to the (small) surviving subspace after prime 1.
//
// Run: magma ro2_search3.m   (long, background)

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf, BadPrimesFromInvariants;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K); F11 := GF(11);
LOGF := "ro2_search3_out.txt"; PrintFile(LOGF, "# ro2_search3" : Overwrite := true);
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

load "mod11_trnm.txt";
prs := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    Append(~prs, <PrimeFromTag(row[2],row[3]), F11!row[4], F11!row[5], row[1]>);
end for;
Sort(~prs, func<a,b|a[4]-b[4]>);

// Part 1: t_2 from g at level Pbad
M0 := HilbertCuspForms(K, Pbad, [2,2]); d0 := Dimension(M0); Id0 := IdentityMatrix(F11,d0);
Vp := VectorSpace(F11,d0); T0 := AssociativeArray();
for pr in prs do
    if Norm(pr[1]+Pbad) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M0, pr[1])), F11); T0[pr[4]] := T;
    chi := ChiG(pr[1])*ChiM1(pr[1]); Vp := Vp meet Kernel(T*T - chi*pr[2]*T + pr[3]*Id0);
    if Dimension(Vp) le 1 and #Keys(T0) ge 24 then break; end if;
end for;
v := Basis(Vp)[1]; t2 := AssociativeArray();
for pr in prs do
    if not IsDefined(T0, pr[4]) then continue; end if;
    w := v*T0[pr[4]]; k := 1; while v[k] eq 0 do k+:=1; end while;
    chi := ChiG(pr[1])*ChiM1(pr[1]); t2[pr[4]] := pr[2] - chi*(w[k]/v[k]);
end for;
T0 := 0;   // free
LOG(Sprintf("Part 1 done: V+(g)=%o, t_2 at %o primes", Dimension(Vp), #Keys(t2)));

// Part 2: Pbad*(2)^3
Level := Pbad*p2^3;
tm := Cputime(); M2 := HilbertCuspForms(K, Level, [2,2]); d2 := Dimension(M2); Id2 := IdentityMatrix(F11,d2);
LOG(Sprintf("Pbad*(2)^3: dim %o (build %o s)", d2, Cputime(tm)));
Ts := []; tps := [];
for pr in prs do
    if Norm(pr[1]+Level) ne 1 or not IsDefined(t2, pr[4]) then continue; end if;
    th := Cputime(); T := ChangeRing(Matrix(HeckeOperator(M2, pr[1])), F11);
    Append(~Ts, T); Append(~tps, pr);
    LOG(Sprintf("  Hecke N=%o done (%o s), stored %o", pr[4], Cputime(th), #Ts));
    if #Ts ge 5 then break; end if;
end for;

LOG("sweeping 8 twists for rho_2 (restriction-optimised)...");
for fi in [1..#funcs] do
    fnc := funcs[fi]; Vp2 := VectorSpace(F11,d2); dead := false;
    for j in [1..#Ts] do
        pr := tps[j]; ev := ChiG(pr[1])*ChiD(pr[1],fnc)*t2[pr[4]];
        if Dimension(Vp2) eq d2 then
            Vp2 := Kernel(Ts[j] - ev*Id2);
        else
            B := BasisMatrix(Vp2); r := Nrows(B);
            Cm := Solution(B, B*Ts[j]) - ev*IdentityMatrix(F11, r);
            kerC := Kernel(Cm);
            if Dimension(kerC) eq 0 then dead := true; break; end if;
            Vp2 := sub< VectorSpace(F11,d2) | [ b*B : b in Basis(kerC) ] >;
        end if;
        if Dimension(Vp2) eq 0 then dead := true; break; end if;
    end for;
    if not dead and Dimension(Vp2) gt 0 then
        LOG(Sprintf("  *** twist #%o: V+ = %o -- rho_2's form at Pbad*(2)^3! ***", fi, Dimension(Vp2)));
    end if;
end for;
LOG("(no surviving twist for rho_2 at Pbad*(2)^3 unless noted above)");
LOG("DONE."); exit;
