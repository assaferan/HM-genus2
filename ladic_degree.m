// ladic_degree.m -- cheap lower bound on the Hecke-field degree [E:Q] of the newform
// matched by the kernel survivor, WITHOUT isolating the eigenform (no NewformDecomposition).
//
// The survivor v is a mod-l eigenvector (v*T_P = t_P*v over F_l, t_P = -#C(F_P) mod l).
// Hensel-lift (v, {a_P}) over Z/l^k with the LEFT eigenvector equation v*T_P = a_P*v.
// a_P(f) lives in E_lambda (lambda|l the fingerprint prime, residue degree 1), and
//   [E:Q] = sum_{lambda|l} [E_lambda:Q_l] >= [E_lambda:Q_l].
// Two outcomes, both certify d = [E:Q] > 1 (i.e. NOT from an elliptic curve):
//   (a) OBSTRUCTION: a_P not in Z_l (E_lambda != Q_l) -> the Z_l lift fails at low precision.
//   (b) RECOGNITION: a_P in Z_l -> recover its min poly via PowerRelation; degree > 1.
//
// Usage:  magma idx:=3 ladic_degree.m           (sweeps e for the 1-dim survivor)
//         magma idx:=1 EEXP:=3 PREC:=60 ladic_degree.m
SetColumns(0);
load "torsion_data.m";
idx := StringToInteger(Sprintf("%o", idx));
if not assigned EEXP then EEXP := -1; end if; EEXP := StringToInteger(Sprintf("%o", EEXP));
if not assigned PREC then PREC := 40; end if; PREC := StringToInteger(Sprintf("%o", PREC));
if not assigned NCOND then NCOND := 14; end if; NCOND := StringToInteger(Sprintf("%o", NCOND));
PBOUND := 120;

row := torsion_data[idx];
d := row[1]; l := row[2]; lab := row[3]; fc := row[4]; hc := row[5]; N := row[6]; cc := row[7];
K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K); Fl := GF(l);
fpoly := R![K| c[1]+c[2]*a : c in fc]; hh := R![K| c[1]+c[2]*a : c in hc];
C := HyperellipticCurve(fpoly, hh);
p2 := Factorization(2*OK)[1][1];
if d eq 2 then base := ideal<OK | cc[1]+cc[2]*a>;
else Cnd := Conductor(C); base := 1*OK;
  for pe in Factorization(Cnd) do if Norm(pe[1]) mod 2 ne 0 then base := base*pe[1]^pe[2]; end if; end for;
end if;

// fingerprint with INTEGER eigenvalues -#C(F_P)
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
  Tint := [ ChangeRing(Matrix(HeckeOperator(M, t[1])), Integers()) : t in fps ];
  Ifl := IdentityMatrix(Fl, dm); cur := VectorSpace(Fl, dm);
  for i in [1..#Tint] do cur := cur meet Kernel(ChangeRing(Tint[i], Fl) - (Fl!(fps[i][2]))*Ifl); end for;
  sd := Dimension(cur);
  printf "  e=%o dim=%o survivor=%o\n", e, dm, sd;
  if sd ne 1 then continue; end if;

  // ---- Hensel lift of (v, {a_i}) over Z/l^PREC, left eigenvector v*T_i = a_i v ----
  v0 := Basis(cur)[1];
  piv := Rep({j : j in [1..dm] | v0[j] ne 0});
  free := [ j : j in [1..dm] | j ne piv ]; nfree := #free; S := #Tint;
  Zk := Integers(l^PREC);
  v := Vector(Zk, [ Integers()!(v0[j]) : j in [1..dm] ]); v := v*(v[piv]^-1);
  avec := [ Zk!(Integers()!(Fl!(fps[i][2]))) : i in [1..S] ];
  TZ := [ ChangeRing(Tint[i], Zk) : i in [1..S] ];
  obstructed := false; mstop := 0;
  for m in [1..PREC-1] do
    Arows := []; w := [];
    for i in [1..S] do
      Mi := TZ[i] - avec[i]*IdentityMatrix(Zk, dm);
      ri := v*TZ[i] - avec[i]*v;
      for c in [1..dm] do
        assert (Integers()!ri[c]) mod (l^m) eq 0;
        coef := [Fl| ];
        for j in free do Append(~coef, Fl!(Integers()!(Mi[j][c]))); end for;
        for ii in [1..S] do Append(~coef, ii eq i select Fl!(-(Integers()!(v[c]))) else Fl!0); end for;
        Append(~Arows, coef);
        Append(~w, Fl!( -((Integers()!ri[c]) div (l^m)) ));
      end for;
    end for;
    Amat := Matrix(Fl, S*dm, nfree+S, Arows);
    cons, u := IsConsistent(Transpose(Amat), Vector(Fl, w));
    if not cons then obstructed := true; mstop := m; break; end if;
    du := Eltseq(u);
    for t in [1..nfree] do v[free[t]] := v[free[t]] + (l^m)*(Zk!(Integers()!(du[t]))); end for;
    for i in [1..S] do avec[i] := avec[i] + (l^m)*(Zk!(Integers()!(du[nfree+i]))); end for;
    mstop := m+1;
  end for;

  if obstructed then
    printf "  RESULT: lift OBSTRUCTED at l-adic digit m=%o  =>  a_P not in Z_l  =>  [E_lambda:Q_l]>1  =>  [E:Q]>1.\n", mstop;
    printf "  ==> the matched newform is NOT rational; sigma does not arise from an elliptic curve. (cheap branch)\n";
  else
    printf "  lift converged to l^%o; recovering min polys (LLL reconstruction):\n", PREC;
    Zx<X> := PolynomialRing(Integers());
    modp := l^PREC;
    // minimal-degree integer relation with a_P as an l-adic root (a_P in [0,l^PREC))
    // find candidate at HALF precision, then VERIFY it vanishes to FULL precision
    // (a spurious lattice relation with huge coefficients fails the full-precision check)
    pf := PREC div 2; mf := l^pf; mfull := modp;
    Recog := function(ap)
      for dd in [1..Min(dm,25)] do
        n := dd+1;
        B := ZeroMatrix(Integers(), n+1, n+1);
        for j in [0..dd] do B[j+1][j+1] := 1; B[j+1][n+1] := mf*(Modexp(ap,j,mf)); end for;
        B[n+1][n+1] := mf*mf;
        L := LLL(B);
        for r in [1..n+1] do
          c := [ L[r][j+1] : j in [0..dd] ];
          if &and[ ci eq 0 : ci in c ] then continue; end if;
          poly := &+[ c[j+1]*X^j : j in [0..dd] ];
          if Degree(poly) lt 1 then continue; end if;
          if (&+[ c[j+1]*Modexp(ap,j,mfull) : j in [0..dd] ]) mod mfull eq 0 then
            return poly div Content(poly);   // verified to full precision -> genuine
          end if;
        end for;
      end for;
      return Zx!0;
    end function;
    maxdeg := 1;
    for i in [1..S] do
      rel := Recog(Integers()!avec[i]);
      printf "    a_{N(P)=%o}: minpoly degree %o%o\n", Norm(fps[i][1]), Degree(rel),
        (Degree(rel) ge 1 and IsIrreducible(rel)) select " (irreducible)" else "";
      if Degree(rel) gt maxdeg then maxdeg := Degree(rel); end if;
    end for;
    printf "  RESULT: max recovered [Q(a_P):Q] = %o  =>  [E:Q] >= %o.  %o\n",
      maxdeg, maxdeg, maxdeg gt 1 select "NOT from an elliptic curve." else "(inconclusive: try more primes/precision)";
  end if;
  exit;
end for;
printf "no 1-dim survivor found in e-range\n";
exit;
