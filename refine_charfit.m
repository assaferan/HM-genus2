// refine_charfit.m -- read saved sign data, test if eps is a quadratic Hecke character.
// Reads refine_serre5_signs.txt (lines: "N p rid  a=.. tTr=.. tNm=..  eps=+/-/0").
// Run: magma refine_charfit.m

SetColumns(0);
K<w> := QuadraticField(5); OK := Integers(K);
QI_norm := 66179;
// rebuild Pbad as the unique prime of norm 66179 (66179 splits in Q(sqrt5))
Pbad := [P[1] : P in Factorization(66179*OK) | Norm(P[1]) eq 66179][1];
p2 := Factorization(2*OK)[1][1];
P3 := Factorization(3*OK)[1][1];
P5 := Factorization(5*OK)[1][1];
f11 := Factorization(11*OK); P11a := f11[1][1]; P11b := f11[2][1];
function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;

LOGF := "refine_charfit_out.txt";
PrintFile(LOGF, "# refine_charfit: is eps a quadratic character?" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

lines := Split(Read("refine_serre5_signs.txt"), "\n");
sd := [];
for ln in lines do
    if #ln eq 0 or ln[1] eq "#" then continue; end if;
    f := [t : t in Split(ln, " ") | #t gt 0];
    p := StringToInteger(f[2]); rid := StringToInteger(f[3]);
    e := f[#f];  // "eps=+" etc
    bit := e[#e];  // last char
    if bit eq "+" then Append(~sd, <PrimeFromTag(p,rid), GF(2)!0>);
    elif bit eq "-" then Append(~sd, <PrimeFromTag(p,rid), GF(2)!1>);
    end if;  // skip ambiguous "0"
end for;
LOG(Sprintf("loaded %o unambiguous sign points", #sd));

supports := [
    <p2^4, "p2^4">,
    <p2^4*P3, "p2^4*3">,
    <p2^4*P5, "p2^4*P5">,
    <p2^4*P11a*P11b, "p2^4*11">,
    <p2^4*Pbad, "p2^4*Pbad">,
    <p2^4*P3*P5, "p2^4*3*P5">,
    <p2^4*P3*P5*P11a*P11b, "p2^4*3*P5*11">,
    <p2^4*P3*P5*P11a*P11b*Pbad, "p2^4*3*P5*11*Pbad">,
    <p2^6*P3*P5*P11a*P11b*Pbad, "p2^6*3*P5*11*Pbad">
];
found := false;
for sp in supports do
  for inf in [[1,2], [Integers()|]] do
    m := sp[1];
    if #inf eq 0 then R, mp := RayClassGroup(m); else R, mp := RayClassGroup(m, inf); end if;
    Q, qm := quo<R | [2*R.j : j in [1..Ngens(R)]]>;
    rk := Ngens(Q);
    if rk eq 0 then continue; end if;
    rows := [[GF(2)!c : c in Eltseq(qm(p[1] @@ mp))] : p in sd];
    bb := Vector(GF(2), [p[2] : p in sd]);
    A := Matrix(GF(2), rows);
    cons, sol := IsConsistent(Transpose(A), bb);
    if cons then
        nontriv := exists{c : c in Eltseq(sol) | c ne GF(2)!0};
        LOG(Sprintf("  FIT: support %o, oo=%o, rk(R/2R)=%o -> %o (functional %o)",
            sp[2], inf, rk, nontriv select "NONTRIVIAL char" else "trivial", Eltseq(sol)));
        found := true;
    else
        LOG(Sprintf("  no: support %o, oo=%o, rk=%o inconsistent", sp[2], inf, rk));
    end if;
  end for;
end for;
if not found then LOG("NO quadratic character of tried support reproduces eps."); end if;

// Independent necessary test: a true quadratic character is multiplicative, so for the
// SPLIT primes P, Pbar over the same rational p, eps(P)*eps(Pbar) = eps((p)) depends only
// on p; and more basically eps should be consistent. Report the raw +/- balance.
nplus := #[x : x in sd | x[2] eq GF(2)!0];
nminus := #[x : x in sd | x[2] eq GF(2)!1];
LOG(Sprintf("raw sign balance: + = %o, - = %o", nplus, nminus));
LOG("DONE.");
exit;
