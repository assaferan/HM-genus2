// refine_serre3.m -- fast, fully-logged diagnostic of the matched space V.
//
// For each test prime P we restrict T_P to V (2x2 matrix C_P over F11) and print:
//   - eigenvalues of C_P (over F11; "irred/F121" if char poly is irreducible),
//   - the target root set {t1,t2} = roots of q_P^+ = y^2 - tTr y + tNm,
//   - whether C_P's eigenvalues lie in {+-t1,+-t2} (the twist-allowed set),
//   - charpoly(C_P) vs q_P^+ / q_P^-  (the single-eigensystem twist test).
// This settles whether V genuinely realizes rho_1 (+) rho_2 (up to a quadratic twist).
//
// Run: magma refine_serre3.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
Rt<y> := PolynomialRing(F11);
load "mod11_trnm.txt";

LOGF := "refine_serre3_out.txt";
PrintFile(LOGF, "# refine_serre3: per-prime eigenvalue diagnostic of V" : Overwrite := true);
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

targets := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    Append(~targets, <PrimeFromTag(row[2],row[3]), F11!row[4], F11!row[5], row[1], row[2]>);
end for;
Sort(~targets, func<a,b | a[4]-b[4]>);

procedure Diagnose(Level, NFIND, NTEST, tag)
    LOG(Sprintf("==== level %o (norm %o) ====", tag, Norm(Level)));
    M := HilbertCuspForms(K, Level, [2,2]);
    d := Dimension(M); Id := IdentityMatrix(F11, d);
    LOG(Sprintf("  cuspidal dim = %o", d));
    cur := VectorSpace(F11, d); lastdim := d; stable := 0; nf := 0;
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        T2 := T*T;
        cur := cur meet (Kernel(T2 - tg[2]*T + tg[3]*Id) + Kernel(T2 + tg[2]*T + tg[3]*Id));
        T := 0; T2 := 0; nf +:= 1;
        if Dimension(cur) eq lastdim then stable +:= 1; else stable := 0; end if;
        lastdim := Dimension(cur);
        LOG(Sprintf("    [findV] N=%o dim=%o", tg[4], Dimension(cur)));
        if Dimension(cur) eq 0 then break; end if;
        if (Dimension(cur) le 3 and stable ge 4) or nf ge NFIND then break; end if;
    end for;
    dV := Dimension(cur);
    LOG(Sprintf("  dim V = %o", dV));
    if dV eq 0 then LOG("  -> no survivor."); return; end if;
    V := cur; BM := BasisMatrix(V);

    nmatch := 0; ntest := 0;
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        C := Solution(BM, BM*T); T := 0;
        cp := CharacteristicPolynomial(C);          // over F11
        // eigenvalues of C
        evC := [r[1] : r in Roots(cp)];
        // target roots of q^+
        qp := y^2 - tg[2]*y + tg[3];
        qm := y^2 + tg[2]*y + tg[3];
        troots := [r[1] : r in Roots(qp)];
        allowed := {r[1] : r in Roots(qp)} join {r[1] : r in Roots(qm)};  // {+-t1,+-t2}
        evIn := forall{e : e in evC | e in allowed} and (#evC eq Degree(cp));
        // single-eigensystem twist test: charpoly = q^+ or q^- ?
        cpMatch := (cp eq qp) select "q+" else ((cp eq qm) select "q-" else "NO");
        if cpMatch ne "NO" then nmatch +:= 1; end if;
        LOG(Sprintf("  N=%o p=%o: charpoly=%o  Cevs=%o (in+-roots:%o)  q+roots=%o tTr=%o tNm=%o  cp-vs-q:%o",
            tg[4], tg[5], cp, evC, evIn, troots, tg[2], tg[3], cpMatch));
        ntest +:= 1;
        if ntest ge NTEST then break; end if;
    end for;
    LOG(Sprintf("  SUMMARY %o: charpoly(C_P) in {q+,q-} at %o/%o primes", tag, nmatch, ntest));
end procedure;

Diagnose(Pbad*p2, 14, 25, "Pbad*p2");
Diagnose(Pbad,    10, 12, "Pbad");
LOG("DONE.");
exit;
