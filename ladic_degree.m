// ladic_degree.m -- lower bound on the Hecke-field degree [E:Q] of the newform matched
// by the kernel survivor, WITHOUT isolating the eigenform (no NewformDecomposition).
//
// Key fact (multiplicity one): at the maximal ideal m cut out by the fingerprint system
// T_P |-> t_P (mod l), the m-adic Hecke module is free of rank 1, so the mod-l GENERALIZED
// eigenspace G has  dim_{F_l} G = [E_lambda : Q_l]  (lambda|l the fingerprint prime, residue
// degree 1).  Since  [E:Q] = sum_{lambda|l} [E_lambda:Q_l] >= [E_lambda:Q_l] = dim G,
//   dim G > 1   =>   [E:Q] > 1   (newform not rational; sigma not from an elliptic curve).
// A 1-dim survivor (exact eigenspace) rules out other congruent forms, so G is pure.
//
// dim G is computed by the increasing chain of NESTED KERNELS
//   E_1 = survivor,   E_{k+1} = { w : w*(T_i - c_i) in E_k for all i },   G = lim E_k,
// i.e. the same F_l linear algebra as the survivor -- it SCALES to the giant dimensions.
// (For dim G = 1, lambda is split/residue-degree-1 and a_P in Z_l; then an l-adic min-poly
// recognition is attempted, but only at small dim -- that branch is O(dim^3) and does NOT scale.)
//
// Usage:  magma idx:=3 ladic_degree.m
//         magma idx:=1 EEXP:=3 ladic_degree.m        (NCOND:=6, RECOG_MAXDIM:=2000, PREC:=400)
SetColumns(0);
load "torsion_data.m";
idx := StringToInteger(Sprintf("%o", idx));
if not assigned EEXP then EEXP := -1; end if; EEXP := StringToInteger(Sprintf("%o", EEXP));
if not assigned NCOND then NCOND := 6; end if; NCOND := StringToInteger(Sprintf("%o", NCOND));
if not assigned PREC then PREC := 400; end if; PREC := StringToInteger(Sprintf("%o", PREC));
if not assigned RECOG_MAXDIM then RECOG_MAXDIM := 2000; end if; RECOG_MAXDIM := StringToInteger(Sprintf("%o", RECOG_MAXDIM));
PBOUND := 120;

row := torsion_data[idx];
d := row[1]; l := row[2]; lab := row[3]; fc := row[4]; hc := row[5]; N := row[6]; cc := row[7];
K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K); Fl := GF(l);
fpoly := R![K| c[1]+c[2]*a : c in fc]; hh := R![K| c[1]+c[2]*a : c in hc];
C := HyperellipticCurve(fpoly, hh); p2 := Factorization(2*OK)[1][1];
if d eq 2 then base := ideal<OK | cc[1]+cc[2]*a>;
else Cnd := Conductor(C); base := 1*OK;
  for pe in Factorization(Cnd) do if Norm(pe[1]) mod 2 ne 0 then base := base*pe[1]^pe[2]; end if; end for; end if;

fps := [];
for pp in PrimesUpTo(PBOUND) do
  if #fps ge NCOND then break; end if;
  if pp eq 2 or pp eq l then continue; end if;
  for tup in Factorization(pp*OK) do
    P := tup[1]; if Norm(base+P) ne 1 then continue; end if;
    kk,red := ResidueClassField(P); Rk := PolynomialRing(kk); ok:=true; Ck:=0;
    try Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(fpoly)], Rk![red(c):c in Coefficients(hh)]);
    catch e; ok:=false; end try;
    if ok then Append(~fps, <P, -#Ck>); end if;
    if #fps ge NCOND then break; end if;
  end for;
end for;
printf "idx %o (%o) l=%o; fingerprint norms %o\n", idx, lab, l, [Norm(t[1]) : t in fps];

elist := EEXP ge 0 select [EEXP] else [0..Valuation(Conductor(C), p2)];
for e in elist do
  Nlev := base * p2^e;
  M := HilbertCuspForms(K, Nlev, [2,2]); dm := Dimension(M);
  if dm eq 0 then continue; end if;
  // HeckeOperator returns a matrix over Q which is NOT always integral (it happens to be for
  // 881/14303, which is why an intermediate ChangeRing(-, Integers()) went unnoticed here; it
  // fails outright on e.g. 24889 and the Q(sqrt 3) levels). Keep the rational matrix and coerce
  // straight to the target ring, exactly as kernel_torsion.m / grh_kernel.m do. Denominators are
  // prime to l, so Q -> F_l and Q -> Z/l^PREC are both well defined; if one ever were not, the
  // coercion errors loudly rather than silently giving a wrong answer.
  TQ := [ Matrix(HeckeOperator(M, t[1])) : t in fps ];
  Ifl := IdentityMatrix(Fl, dm); S := #fps;
  Mmats := [ ChangeRing(TQ[i], Fl) - (Fl!(fps[i][2]))*Ifl : i in [1..S] ];
  Vs := VectorSpace(Fl, dm);
  E1 := Vs; for i in [1..S] do E1 := E1 meet Kernel(Mmats[i]); end for;
  sd := Dimension(E1);
  printf "  e=%o dim=%o survivor=%o\n", e, dm, sd;
  if sd ne 1 then continue; end if;

  // ---- generalized eigenspace dimension via nested kernels (scalable) ----
  Ecur := E1;
  printf "  building generalized eigenspace (nested kernels):\n";
  repeat
    dprev := Dimension(Ecur);
    Bk := BasisMatrix(Ecur);
    Enext := Vs;
    for i in [1..S] do
      Kst := Kernel(VerticalJoin(Mmats[i], -Bk));   // {(w,mu): w*M_i - mu*Bk = 0} = {w: w*M_i in E_k}
      prj := sub< Vs | [ Vs ! [ b[j] : j in [1..dm] ] : b in Basis(Kst) ] >;
      Enext := Enext meet prj;
    end for;
    Ecur := Enext;
    printf "    dim = %o\n", Dimension(Ecur);
  until Dimension(Ecur) eq dprev;
  dimG := Dimension(Ecur);
  printf "  GENERALIZED-EIGENSPACE dim = [E_lambda:Q_l] = %o\n", dimG;

  if dimG gt 1 then
    printf "  RESULT: [E:Q] >= %o > 1  =>  matched newform NOT rational; sigma does not arise from an elliptic curve.\n", dimG;
    printf "LADIC-DEG %-11o idx=%o l=%o e=%o dim=%o dimG=%o degree>=%o NOT-ELLIPTIC\n", lab, idx, l, e, dm, dimG, dimG;
    exit;
  end if;

  printf "  dim G = 1 (lambda split, residue degree 1): a_P in Z_l -- generalized-eigenspace test inconclusive.\n";
  if dm gt RECOG_MAXDIM then
    printf "  RESULT: INCONCLUSIVE at dim %o (l-adic recognition is O(dim^3), does not scale; would need eigenform isolation).\n", dm;
    printf "LADIC-DEG %-11o idx=%o l=%o e=%o dim=%o dimG=1 degree>=1 INCONCLUSIVE-SPLIT-BIGDIM\n", lab, idx, l, e, dm;
    exit;
  end if;

  // ---- small-dim fallback: l-adic Hensel lift + LLL recognition of a_P's min poly ----
  printf "  dim %o <= %o: l-adic min-poly recognition (PREC=%o)...\n", dm, RECOG_MAXDIM, PREC;
  v0 := Basis(E1)[1]; piv := Rep({j : j in [1..dm] | v0[j] ne 0});
  free := [ j : j in [1..dm] | j ne piv ]; nfree := #free;
  Zk := Integers(l^PREC);
  v := Vector(Zk, [ Integers()!(v0[j]) : j in [1..dm] ]); v := v*(v[piv]^-1);
  avec := [ Zk!(Integers()!(Fl!(fps[i][2]))) : i in [1..S] ];
  TZ := [ ChangeRing(TQ[i], Zk) : i in [1..S] ];
  for m in [1..PREC-1] do
    Arows := []; w := [];
    for i in [1..S] do
      Mi := TZ[i] - avec[i]*IdentityMatrix(Zk, dm); ri := v*TZ[i] - avec[i]*v;
      for c in [1..dm] do
        coef := [Fl| ]; for j in free do Append(~coef, Fl!(Integers()!(Mi[j][c]))); end for;
        for ii in [1..S] do Append(~coef, ii eq i select Fl!(-(Integers()!(v[c]))) else Fl!0); end for;
        Append(~Arows, coef); Append(~w, Fl!( -((Integers()!ri[c]) div (l^m)) ));
      end for;
    end for;
    cons, u := IsConsistent(Transpose(Matrix(Fl, S*dm, nfree+S, Arows)), Vector(Fl, w));
    assert cons; du := Eltseq(u);
    for t in [1..nfree] do v[free[t]] := v[free[t]] + (l^m)*(Zk!(Integers()!(du[t]))); end for;
    for i in [1..S] do avec[i] := avec[i] + (l^m)*(Zk!(Integers()!(du[nfree+i]))); end for;
  end for;
  Zx<X> := PolynomialRing(Integers()); modp := l^PREC; mf := l^(PREC div 2);
  Recog := function(ap)
    for dd in [1..25] do
      n := dd+1; B := ZeroMatrix(Integers(), n+1, n+1);
      for j in [0..dd] do B[j+1][j+1]:=1; B[j+1][n+1]:=mf*(Modexp(ap,j,mf)); end for;
      B[n+1][n+1]:=mf*mf; L := LLL(B);
      for r in [1..n+1] do
        c := [L[r][j+1] : j in [0..dd]];
        if &and[ci eq 0: ci in c] then continue; end if;
        poly := &+[c[j+1]*X^j : j in [0..dd]];
        if Degree(poly) ge 1 and (&+[c[j+1]*Modexp(ap,j,modp):j in [0..dd]]) mod modp eq 0 then return poly div Content(poly); end if;
      end for;
    end for; return Zx!0;
  end function;
  maxdeg := 1;
  for i in [1..S] do
    rec := Recog(Integers()!avec[i]);
    if Degree(rec) gt maxdeg then maxdeg := Degree(rec); end if;
  end for;
  printf "  RESULT: max recovered [Q(a_P):Q] = %o  =>  [E:Q] >= %o.  %o\n",
    maxdeg, maxdeg, maxdeg gt 1 select "NOT from an elliptic curve." else "(inconclusive: more primes/precision)";
  printf "LADIC-DEG %-11o idx=%o l=%o e=%o dim=%o dimG=1 degree>=%o %o\n",
    lab, idx, l, e, dm, maxdeg, maxdeg gt 1 select "NOT-ELLIPTIC-VIA-RECOG" else "INCONCLUSIVE-RECOG";
  exit;
end for;
printf "no 1-dim survivor found in e-range\n";
printf "LADIC-DEG %-11o idx=%o l=%o NO-SURVIVOR\n", lab, idx, l;
exit;
