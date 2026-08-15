// emit_cutters.m -- generator: emit the Magma Hecke data (Hecke field + cutters) for one
// matched form, in a form that loads back into Magma. Reusable for any curve in the dataset.
//
//   magma lab:="14303.1" e:=0 emit_cutters.m         // prints the block to stdout
//   magma lab:="881.1"   e:=3 out:="hecke_cutters.m" emit_cutters.m   // appends to a file
//
// `lab` = a label in torsion_data.m; `e` = the 2-part exponent found by the sweep (level =
// prime-to-2 conductor * p2^e). Optional `out` appends the block to that file; else stdout.
// Optional `pbound` (default 40) bounds the cutter primes.
SetColumns(0);
load "torsion_data.m";
if not assigned pbound then pbound := 40; end if;
if not assigned e then e := 0; end if;
if Type(e) eq MonStgElt then e := StringToInteger(e); end if;            // cmdline args arrive as strings
if Type(pbound) eq MonStgElt then pbound := StringToInteger(pbound); end if;

row := 0;
for r in torsion_data do if r[3] eq lab then row := r; break; end if; end for;
error if row cmpeq 0, "emit_cutters: label not found in torsion_data.m: " cat lab;
d := row[1]; l := row[2]; fc := row[4]; hc := row[5]; cc := row[7];
K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K); Fl := GF(l);
Qy<y> := PolynomialRing(Rationals());
f := R![K| c[1]+c[2]*a : c in fc]; hh := R![K| c[1]+c[2]*a : c in hc]; C := HyperellipticCurve(f, hh);
p2 := Factorization(2*OK)[1][1];
if d eq 2 then base := ideal<OK | cc[1]+cc[2]*a>;
else base := 1*OK; for pe in Factorization(Conductor(C)) do if Norm(pe[1]) mod 2 ne 0 then base := base*pe[1]^pe[2]; end if; end for; end if;
N := base * p2^(Integers()!e);

// fingerprint to select the matching orbit
fp := [];
for pp in PrimesUpTo(pbound) do
    if pp eq 2 or pp eq l then continue; end if;
    for tup in Factorization(pp*OK) do
        P := tup[1]; if Norm(N+P) ne 1 then continue; end if;
        k,red := ResidueClassField(P); Rk := PolynomialRing(k); ok := true; Ck := 0;
        try Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(f)], Rk![red(c):c in Coefficients(hh)]);
        catch ee; ok := false; end try;
        if ok then Append(~fp, <P, Fl!(-#Ck)>); end if;
    end for;
end for;

Mnew := NewSubspace(HilbertCuspForms(K, N, [2,2]));
decomp := NewformDecomposition(Mnew);
matchnf := 0;
for nf in decomp do
    eig := Eigenform(nf); Ff := HeckeEigenvalueField(nf); OF := Integers(Ff);
    for lam in [pr[1]:pr in Factorization(l*OF)] do
        kf,m := ResidueClassField(lam); good := true;
        for pr in fp do if m(OF!HeckeEigenvalue(eig,pr[1])) ne kf!(Integers()!pr[2]) then good := false; break; end if; end for;
        if good then matchnf := nf; break; end if;
    end for;
    if matchnf cmpne 0 then break; end if;
end for;
error if matchnf cmpeq 0, "emit_cutters: no matching orbit at level norm " cat IntegerToString(Norm(N)) cat " for " cat lab;

eig := Eigenform(matchnf); Ff := HeckeEigenvalueField(matchnf);
tag := "";
for ch in Eltseq(lab) do tag := tag cat (ch eq "." select "_" else ch); end for;
// Hecke field invariants (comment line)
Kh := NumberField(DefiningPolynomial(Ff)); dh := Discriminant(MaximalOrder(Kh));
lmfdb := (Degree(Kh) le 6) select Sprintf("try LMFDB nf disc %o deg %o", AbsoluteValue(dh), Degree(Kh)) else "not in LMFDB (large degree/disc)";

lines := [
  Sprintf("// ---- %o  (F=Q(sqrt%o), l=%o; level norm %o = p2^%o * (prime-to-2 cond); orbit dim %o) ----",
          lab, d, l, Norm(N), e, Dimension(matchnf)),
  Sprintf("K<a> := QuadraticField(%o); OK := Integers(K);   // self-contained: a = sqrt(%o)", d, d),
  Sprintf("//      Hecke field: deg %o, disc %o, Galois %oT%o; %o",
          Degree(Kh), dh, Degree(Kh),
          (Degree(Kh) le 15) select TransitiveGroupIdentification(GaloisGroup(DefiningPolynomial(Ff))) else 0, lmfdb),
  Sprintf("heckefield_%o := %o;", tag, Qy!DefiningPolynomial(Ff))
];
cut := Sprintf("cutters_%o := [", tag);
cnt := 0;
for pp in PrimesUpTo(pbound) do
    for tup in Factorization(pp*OK) do
        P := tup[1]; if Norm(N+P) ne 1 then continue; end if;
        isp, g := IsPrincipal(P);
        if isp then idgen := Sprintf("%o", K!g);   // class number 1: single generator
        else gg := [K!x : x in Generators(P)];     // non-principal: emit the full generator list
             idgen := &cat[ Sprintf("%o%o", gg[i], i lt #gg select ", " else "") : i in [1..#gg] ]; end if;
        ap := HeckeEigenvalue(eig, P);
        cut := cut cat Sprintf("\n  < ideal<OK | %o>, %o >,   // N(P) = %o", idgen, Qy!MinimalPolynomial(ap), Norm(P));
        cnt +:= 1;
    end for;
    if cnt ge 6 then break; end if;
end for;
cut := cut cat "\n];";
Append(~lines, cut);
block := "\n" cat &cat[ ln cat "\n" : ln in lines ];

if assigned out then PrintFile(out, block); printf "appended %o to %o\n", lab, out;
else printf "%o", block; end if;
quit;
