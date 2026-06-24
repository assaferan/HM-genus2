// extract_eigensystems.m -- (a) identify the quadratic twist eps and extract the actual
// Hecke eigensystems of the forms realizing rho_1, rho_2; (b) verify over many primes.
//
// V (dim 2, in S_2(p*(2)) (x) F_11) is the eigenform span. For each prime P we restrict the
// Hecke operator T_P to V (a 2x2 matrix X_P) and read its eigenvalues = {a_P(f_1),a_P(f_2)}.
// Compare to the fingerprint roots {t,t'} (roots of y^2 - Tr_P y + Nm_P):
//   {a_P} = {t,t'}  => eps(P)=+1 ;  {a_P} = {-t,-t'} => eps(P)=-1 ;  else mixed/fail.
// Then identify eps among quadratic characters of K=Q(sqrt5) ramified only at 2.
//
// Run: magma extract_eigensystems.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
R<y> := PolynomialRing(F11);
load "mod11_trnm.txt";
plog := "extract_progress.txt"; PrintFile(plog, "# extract" : Overwrite := true);
procedure LOG(s) PrintFile(plog, s); end procedure;

QI_A := [K| 1, 1/1459240*(243125*w - 482787),
    1/2229718720*(-54798934*w + 127710391),
    1/8517525510400*(182422196780*w - 406668692089),
    1/65073894899456000*(38134182372761*w - 85270624049895) ];
Pbad := [P : P in BadPrimesFromInvariants(QI_A, K) | Norm(P) eq 66179][1];
p2 := Factorization(2*OK)[1][1];
Level := Pbad*p2;
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

printf "Building space at level norm %o...\n", Norm(Level);
M := HilbertCuspForms(K, Level, [2,2]);
d := Dimension(M); Id := IdentityMatrix(F11, d);
printf "  dim = %o\n", d;

// ---- build V over the first several coprime primes (twist-robust) ----
cur := VectorSpace(F11, d); built := 0;
for tg in targets do
    if Norm(tg[1] + Level) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11); T2 := T*T;
    kp := Kernel(T2 - tg[2]*T + tg[3]*Id); km := Kernel(T2 + tg[2]*T + tg[3]*Id);
    cur := cur meet (kp + km); T := 0; T2 := 0; kp := 0; km := 0;
    built +:= 1; LOG(Sprintf("build prime N=%o: dim V=%o", tg[4], Dimension(cur)));
    if built ge 8 then break; end if;
end for;
printf "  dim V = %o (built from %o primes)\n", Dimension(cur), built;
B := BasisMatrix(cur);   // 2 x d

// ---- extract eigensystems + twist, verifying over many primes ----
printf "\n%-7o %-4o %-12o %-12o %o\n", "N(P)", "p", "{a_P(f_i)}", "{t,t'}", "eps(P)";
printf "%o\n", "-"^58;
epsdata := [];   // <P-ideal, eps>  for character identification
nver := 0; nfail := 0; nmixed := 0;
cnt := 0;
for tg in targets do
    if Norm(tg[1] + Level) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
    BT := B*T;
    fl, X := IsConsistent(B, BT);   // solve X*B = B*T  (V Hecke-stable => consistent)
    if not fl then printf "%-7o %-4o  V not T_P-stable!\n", tg[4], tg[5]; nfail +:= 1; T:=0; continue; end if;
    evsX := { e[1] : e in Eigenvalues(X) };
    qroots := { r[1] : r in Roots(y^2 - tg[2]*y + tg[3]) };
    negq := { -r : r in qroots };
    inpos := evsX subset qroots; inneg := evsX subset negq;
    if inpos and inneg then epstr := "amb"; nver +:= 1;       // a_P in both (e.g. 0 / symmetric)
    elif inpos then epstr := "+1"; nver +:= 1; Append(~epsdata, <tg[1], 1, tg[5], tg[6]>);
    elif inneg then epstr := "-1"; nver +:= 1; Append(~epsdata, <tg[1], -1, tg[5], tg[6]>);
    else epstr := "FAIL"; nfail +:= 1; end if;
    printf "%-7o %-4o %-10o %-12o eps=%o\n", tg[4], tg[5], evsX, qroots, epstr;
    T := 0;
    cnt +:= 1; if cnt ge 24 then break; end if;
end for;
printf "\nverified (a_P in +-{q-roots}): %o ; FAIL: %o\n", nver, nfail;

// ---- identify eps: a quadratic character of K ramified only at 2 ----
// candidates K(sqrt(d)): test eps_d(P) = +1 if P splits in K(sqrt d), else -1.
printf "\nIdentifying eps among quadratic chars of K ramified at 2:\n";
// save (p, rid, eps) for offline use
PrintFile("eps_data.txt", "// <p, rid, eps(P)> from extract_eigensystems.m" : Overwrite := true);
PrintFile("eps_data.txt", Sprintf("epsraw := %o;", [<e[3],e[4],e[2]> : e in epsdata]));
printf "saved %o (P,eps) points to eps_data.txt\n", #epsdata;

// identify eps as an order-2 Hecke character: sweep moduli, report best match + conductor.
p2 := Factorization(2*OK)[1][1];
function BestChar(m, inf)
    Rcl, mp := RayClassGroup(m, inf); ng := Ngens(Rcl);
    ords := [ Order(Rcl.i) : i in [1..ng] ];
    choices := [ (IsEven(ords[i]) select [1,-1] else [1]) : i in [1..ng] ];
    best := 0; bestimg := [];
    for tup in CartesianProduct(choices) do
        img := [ tup[i] : i in [1..ng] ]; c := 0;
        for e in epsdata do
            es := Eltseq(e[1] @@ mp);
            if &*[ img[i]^es[i] : i in [1..ng] ] eq e[2] then c +:= 1; end if;
        end for;
        if c gt best then best := c; bestimg := img; end if;
    end for;
    return best, Rcl, mp, bestimg, ng;
end function;
printf "\nIdentifying eps (order-2 Hecke char), %o data points:\n", #epsdata;
for ms in [<"(2)^4(5)(66179)", p2^4*(5*OK)*(66179*OK)>, <"(2)^6(5)(66179)", p2^6*(5*OK)*(66179*OK)>,
           <"(2)^4(5)", p2^4*(5*OK)>, <"(2)^4(66179)", p2^4*(66179*OK)>] do
    for inf in [[1,2],[Integers()|]] do
        b := BestChar(ms[2], inf);
        printf "  modulus %-16o inf=%o: best matches %o/%o\n", ms[1], inf, b, #epsdata;
    end for;
end for;
