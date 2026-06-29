// refine_serre5.m -- pin down (or refute) the twist character for the Pbad eigenform.
//
// At level Pbad there is a unique eigenform f with a_P(f) in {+-t1,+-t2} at every prime.
// Define eps(P) = +1 if a_P in {t1,t2} (root of q^+), -1 if in {-t1,-t2} (root of q^-).
// f witnesses "rho_i modular up to a quadratic twist" iff eps is a genuine quadratic Hecke
// character.  We (a) dump the raw <N,p,a_P,t1,t2,eps> for ~120 primes, and (b) test whether
// eps is a quadratic character of conductor supported in {p2,3,5,11,66179, oo} (broader than
// before).  If it fits -> identify conductor; if not -> the relation is NOT a quadratic twist.
//
// Run: magma refine_serre5.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
load "mod11_trnm.txt";

LOGF := "refine_serre5_out.txt";
PrintFile(LOGF, "# refine_serre5: twist character analysis" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;
DATF := "refine_serre5_signs.txt";
PrintFile(DATF, "# N p rid a_P t1 t2 eps" : Overwrite := true);

QI_A := [K| 1, 1/1459240*(243125*w - 482787),
    1/2229718720*(-54798934*w + 127710391),
    1/8517525510400*(182422196780*w - 406668692089),
    1/65073894899456000*(38134182372761*w - 85270624049895) ];
Pbad := [P : P in BadPrimesFromInvariants(QI_A, K) | Norm(P) eq 66179][1];
p2 := Factorization(2*OK)[1][1];
function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;

targets := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    Append(~targets, <PrimeFromTag(row[2],row[3]), F11!row[4], F11!row[5], row[1], row[2], row[3]>);
end for;
Sort(~targets, func<a,b | a[4]-b[4]>);

NFIND := 14; NTEST := 120;
M := HilbertCuspForms(K, Pbad, [2,2]);
d := Dimension(M); Id := IdentityMatrix(F11, d);
LOG(Sprintf("level Pbad norm %o, cuspidal dim %o", Norm(Pbad), d));
V := VectorSpace(F11, d); nf := 0;
for tg in targets do
    if Norm(tg[1] + Pbad) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11); T2 := T*T;
    V := V meet (Kernel(T2 - tg[2]*T + tg[3]*Id) + Kernel(T2 + tg[2]*T + tg[3]*Id));
    T := 0; T2 := 0; nf +:= 1;
    if nf ge NFIND then break; end if;
end for;
LOG(Sprintf("dim Vtw = %o", Dimension(V)));
assert Dimension(V) eq 1;
BM := BasisMatrix(V);

// collect signs over NTEST primes
sd := [];           // <ideal, GF2 bit>  (unambiguous primes)
nt := 0; fails := 0;
for tg in targets do
    if Norm(tg[1] + Pbad) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
    C := Solution(BM, BM*T); T := 0;
    a := C[1,1];                 // dim V = 1, so C is 1x1 = the eigenvalue
    inP := (a^2 - tg[2]*a + tg[3]) eq F11!0;
    inM := (a^2 + tg[2]*a + tg[3]) eq F11!0;
    // roots of q+ for display
    epsstr := "?";
    if inP and not inM then epsstr := "+"; Append(~sd, <tg[1], GF(2)!0>);
    elif inM and not inP then epsstr := "-"; Append(~sd, <tg[1], GF(2)!1>);
    elif inP and inM then epsstr := "0"; // ambiguous
    else epsstr := "FAIL"; fails +:= 1; end if;
    PrintFile(DATF, Sprintf("%o %o %o  a=%o tTr=%o tNm=%o  eps=%o",
        tg[4], tg[5], tg[6], a, tg[2], tg[3], epsstr));
    nt +:= 1;
    if nt ge NTEST then break; end if;
end for;
LOG(Sprintf("tested %o primes, FAIL=%o, sign points=%o", nt, fails, #sd));

// broad quadratic-character search
P3 := Factorization(3*OK)[1][1];
P5 := Factorization(5*OK)[1][1];
f11 := Factorization(11*OK); P11a := f11[1][1]; P11b := f11[2][1];
LOG("character search over support {p2^<=4, 3, 5, 11(both), Pbad}, with/without oo:");
found := false;
// try a sequence of moduli of increasing support
supports := [
    <p2^4, "p2^4">,
    <p2^4*P3, "p2^4*3">,
    <p2^4*P5, "p2^4*P5">,
    <p2^4*P11a*P11b, "p2^4*11">,
    <p2^4*Pbad, "p2^4*Pbad">,
    <p2^4*P3*P5*P11a*P11b, "p2^4*3*P5*11">,
    <p2^4*P3*P5*P11a*P11b*Pbad, "p2^4*3*P5*11*Pbad">
];
for sp in supports do
  for inf in [[1,2], []] do
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
        LOG(Sprintf("  FIT: support %o, oo=%o, rk=%o -> %o (func %o)",
            sp[2], inf, rk, nontriv select "NONTRIVIAL char" else "trivial", Eltseq(sol)));
        found := true;
    end if;
  end for;
end for;
if not found then
    LOG("  NO quadratic character (any tried support) reproduces the signs.");
    LOG("  => the eigenform's relation to rho_i is NOT a global quadratic twist.");
end if;
LOG("DONE.");
exit;
