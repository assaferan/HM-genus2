// kernel_w24.m -- kernel method at RAISED weight [2,4], NO decomposition.
// Hecke operators live over the weight field E_wt = Q(zeta_8); reduce mod a prime
// above 5 (residue F_25) and intersect ker(T_(l) - c_l) over inert primes, where
// the induced/normalized inert condition is
//    a_(l)(f) = -a2(A,l) / Norm((l))^d  (mod 5),  d=1 for [2,4],  Norm((l))=l^2.
// Nonzero surviving eigenspace = match.  Discrimination control: c_l -> c_l+1 -> 0.
// Params: idx.  Run: magma -b idx:=20 kernel_w24.m
SetColumns(0);
idx := StringToInteger(idx); k1 := StringToInteger(k1); k2 := StringToInteger(k2);
load "goal1_data.m";
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
ent := [e : e in entries | e[1] eq idx][1];
d := ent[2]; serre := ent[4]; fc := ent[5]; hc := ent[6];
F<w> := QuadraticField(d); OF := Integers(F); F5 := GF(5);
C := HyperellipticCurve(&+[fc[i]*x^(i-1):i in [1..#fc]], &+[hc[i]*x^(i-1):i in [1..#hc]]);
dC := Integers()!Discriminant(C);
J := ideal<OF|1>;
for pe in Factorization(serre) do
    P := Factorization(pe[1]*OF)[1][1]; J := J * P^(pe[2] div Valuation(Norm(P), pe[1]));
end for;
lev := J;
dd := (Max(k1,k2)-2) div 2;

// inert conditions from A: <(l) ideal, c_l in F5, l>
conds := [];
for l in PrimesInInterval(3, 120) do
    if #conds ge 8 then break; end if;
    if l eq 5 or (d mod l eq 0) or (dC mod l eq 0) then continue; end if;
    if IsSquare(GF(l)!d) then continue; end if;          // want inert
    if Norm(lev) mod l eq 0 then continue; end if;
    LA := LPolynomial(ChangeRing(C, GF(l)));
    b := Coefficient(LA,2);                                // a2(A,l)
    // a_(l)(f) = -b / (l^2)^d  mod 5
    c := (F5!(-b)) / ((F5!l)^(2*dd));
    Append(~conds, < ideal<OF|l>, c, l >);
end for;

t0 := Cputime();
M := HilbertCuspForms(F, lev, [k1,k2]); dm := Dimension(M);
printf "idx %o: level norm %o, [%o,%o] dim %o (built %.1os)\n", idx, Norm(lev), k1, k2, dm, Cputime(t0);
if dm eq 0 then printf "NO-MATCH (empty)\n"; exit; end if;
if dm gt 12000 then printf "SKIP-BIG dim %o (kernel over weight field too big)\n", dm; exit; end if;

// weight field and a prime above 5
T0 := HeckeOperator(M, conds[1][1]);
Ewt := BaseRing(Parent(T0)); Owt := Integers(Ewt);
lam := Factorization(5*Owt)[1][1]; Fq, red := ResidueClassField(lam);
printf "weight field deg %o; residue field F_%o\n", Degree(Ewt), #Fq;
function Reduce(T)
    es := Eltseq(T);
    Dl := LCM([ Denominator(e) : e in es ]);   // denominators are powers of 2, coprime to 5
    error if GCD(Dl,5) ne 1, "denominator divisible by 5";
    return (Fq!Dl)^-1 * Matrix(Fq, dm, dm, [ red(Owt!(Dl*e)) : e in es ]);
end function;

V := VectorSpace(Fq, dm); I := IdentityMatrix(Fq, dm);
cur := V; cur2 := V;    // true and control (c+1)
for cc in conds do
    Tp := Reduce(HeckeOperator(M, cc[1]));
    cval := Fq!cc[2];
    cur  := cur  meet Kernel(Tp - cval*I);
    cur2 := cur2 meet Kernel(Tp - (cval+1)*I);
    printf "  l=%o: true cum %o, ctrl cum %o\n", cc[3], Dimension(cur), Dimension(cur2);
    if Dimension(cur) eq 0 and Dimension(cur2) eq 0 then break; end if;
end for;
status := (Dimension(cur) gt 0 and Dimension(cur2) eq 0) select "MATCH" else
          ((Dimension(cur) gt 0) select "MATCH?(ctrl nonzero)" else "NO-MATCH");
printf ">>> idx %o [%o,%o] kernel: %o (true %o / ctrl %o)  [%.1os]\n",
    idx, k1, k2, status, Dimension(cur), Dimension(cur2), Cputime(t0);
exit;
