// refine_serre.m -- refinements of the mod-11 Serre match for Curve A / Q(sqrt5).
//
// Builds on match_hmf_twist.m. Three refinements in one pass:
//  (1) TWIST: the matched space V (dim 2, Hecke-stable) lets us restrict each T_P to a
//      2x2 matrix C_P over F11. Then det C_P = Nm_P (consistency at every prime) and
//      trace C_P = eps(P) * tTr_P, so eps(P) = (trace C_P)/tTr_P in {+1,-1} reads off the
//      quadratic twist prime-by-prime. We then identify eps as an explicit quadratic
//      ray-class character of 2-power conductor.
//  (2) VERIFY OVER MORE PRIMES: do the det/trace check at NMAX primes (>> the original 8).
//  (3) LEVEL EXPONENT: confirm dim V = 0 at level Pbad (no 2-part) and = 2 at Pbad*(2),
//      so the 2-part of the conductor is exactly p2^1.
//
// Run: magma refine_serre.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
load "mod11_trnm.txt";   // trnm := [ <Np,p,rid,Tr,Nm>, ... ]

LOGF := "refine_serre_out.txt";
PrintFile(LOGF, "# refine_serre: twist + verify + level-exponent" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

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

// targets sorted by ascending norm (cheapest Hecke ops first), coprime to 2*11
targets := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    Append(~targets, <PrimeFromTag(row[2],row[3]), F11!row[4], F11!row[5], row[1], row[2]>);
end for;
Sort(~targets, func<a,b | a[4]-b[4]>);

Level := Pbad*p2;
LOG(Sprintf("level = Pbad*p2, norm %o", Norm(Level)));
M := HilbertCuspForms(K, Level, [2,2]);
d := Dimension(M);
LOG(Sprintf("cuspidal dim = %o", d));
Id := IdentityMatrix(F11, d);

// ---- PASS 1: find V (twist-robust intersection) on small primes ----
LOG("PASS 1: locating V ...");
cur := VectorSpace(F11, d); lastdim := d; stable := 0; V := cur;
for tg in targets do
    if Norm(tg[1] + Level) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
    T2 := T*T;
    kp := Kernel(T2 - tg[2]*T + tg[3]*Id);
    km := Kernel(T2 + tg[2]*T + tg[3]*Id);
    cur := cur meet (kp + km);
    T := 0; T2 := 0; kp := 0; km := 0;
    if Dimension(cur) eq lastdim then stable +:= 1; else stable := 0; end if;
    lastdim := Dimension(cur);
    LOG(Sprintf("  N=%o: dim(V)=%o (stable %o)", tg[4], Dimension(cur), stable));
    if Dimension(cur) eq 0 then LOG("  V collapsed -- abort"); exit; end if;
    if Dimension(cur) le 2 and stable ge 4 then V := cur; break; end if;
end for;
LOG(Sprintf("PASS 1 done: dim V = %o", Dimension(V)));
BM := BasisMatrix(V);   // (dim V) x d

// ---- PASS 2: verify det/trace and extract twist sign over NMAX primes ----
NMAX := 32;
LOG(Sprintf("PASS 2: restricting T_P to V over up to %o primes ...", NMAX));
// store <P_ideal, b in GF(2)> for eps identification; b=0 if sign +1, 1 if -1
signdata := [];
nok := 0; nbad := 0; nused := 0;
for tg in targets do
    if Norm(tg[1] + Level) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
    C := Solution(BM, BM*T);      // restriction of T_P to V, (dim V)x(dim V)
    T := 0;
    dt := Determinant(C); tr := Trace(C);
    detok := (dt eq tg[3]);
    if detok then nok +:= 1; else nbad +:= 1; end if;
    sgn := "?";
    if tg[2] ne 0 then
        r := tr / tg[2];                 // expected eps(P) in {1,-1}
        if r eq F11!1 then sgn := "+"; Append(~signdata, <tg[1], GF(2)!0>);
        elif r eq F11!(-1) then sgn := "-"; Append(~signdata, <tg[1], GF(2)!1>);
        else sgn := Sprintf("BAD(%o)", r); end if;
    end if;
    LOG(Sprintf("  N=%o p=%o: det=%o (Nm=%o %o)  tr=%o tTr=%o  eps=%o",
        tg[4], tg[5], dt, tg[3], detok select "OK" else "MISMATCH", tr, tg[2], sgn));
    nused +:= 1;
    if nused ge NMAX then break; end if;
end for;
LOG(Sprintf("PASS 2 done: det matched %o/%o primes (%o mismatch); %o sign points",
    nok, nused, nbad, #signdata));

// ---- identify eps as a quadratic ray-class character of 2-power conductor ----
LOG("IDENTIFY eps: searching quadratic characters of conductor p2^e ...");
found := false;
for e in [1..8] do
    R, mp := RayClassGroup(p2^e, [1,2]);   // include both real places
    // quotient by squares: elementary abelian 2-group Q = R/2R
    Q, qm := quo<R | [2*R.i : i in [1..Ngens(R)]]>;
    rk := Ngens(Q);
    if rk eq 0 then continue; end if;
    // build linear system over GF(2): rows = sign points, cols = generators of Q
    rows := []; rhs := [];
    for sd in signdata do
        x := qm(sd[1] @@ mp);
        Append(~rows, [GF(2)!c : c in Eltseq(x)]);
        Append(~rhs, sd[2]);
    end for;
    A := Matrix(GF(2), rows);            // (#pts) x rk
    b := Vector(GF(2), rhs);
    consistent, sol := IsConsistent(Transpose(A), b);  // f*A^T = b  => A*f^T... solve f with A*f = b
    if consistent then
        // verify the recovered character is NONTRIVIAL (otherwise eps trivial = no twist)
        nontriv := exists{c : c in Eltseq(sol) | c ne GF(2)!0};
        LOG(Sprintf("  e=%o: rk(R/2R)=%o  CONSISTENT, eps %o (functional %o)",
            e, rk, nontriv select "NONTRIVIAL" else "TRIVIAL", Eltseq(sol)));
        if nontriv and not found then
            found := true;
            LOG(Sprintf("  *** eps identified: quadratic char of conductor dividing p2^%o ***", e));
            // describe: which class-group generators it is nontrivial on
        end if;
    else
        LOG(Sprintf("  e=%o: rk(R/2R)=%o  inconsistent", e, rk));
    end if;
end for;
if not found then LOG("eps: no nontrivial quadratic char of 2-power conductor fits (check signs)"); end if;

// ---- LEVEL EXPONENT: V at Pbad (no 2-part) should be 0 ----
LOG("LEVEL EXPONENT: checking dim V at level Pbad (norm 66179) ...");
M0 := HilbertCuspForms(K, Pbad, [2,2]);
d0 := Dimension(M0);
Id0 := IdentityMatrix(F11, d0);
cur0 := VectorSpace(F11, d0); used0 := 0;
for tg in targets do
    if Norm(tg[1] + Pbad) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M0, tg[1])), F11);
    T2 := T*T;
    kp := Kernel(T2 - tg[2]*T + tg[3]*Id0);
    km := Kernel(T2 + tg[2]*T + tg[3]*Id0);
    cur0 := cur0 meet (kp + km);
    T := 0; T2 := 0; kp := 0; km := 0;
    used0 +:= 1;
    LOG(Sprintf("  [Pbad] N=%o: dim=%o", tg[4], Dimension(cur0)));
    if Dimension(cur0) eq 0 then break; end if;
    if used0 ge 8 then break; end if;
end for;
LOG(Sprintf("LEVEL EXPONENT: dim V at Pbad = %o (0 => 2-part exponent is exactly 1)", Dimension(cur0)));

LOG("DONE.");
exit;
