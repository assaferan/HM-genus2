// refine_serre4.m -- the CORRECT test: single eigenform realizing rho_i up to a quadratic twist.
//
// A Hilbert newform f matches rho_i iff a_P(f) = tr rho_i(Frob_P) = a root of q_P^+ at every P
// (untwisted), or a_P(f) = eps(P)*root, i.e. a root of q_P^+ (eps=+1) or q_P^- (eps=-1) with
// eps a fixed quadratic character (twist).  So the candidate is a 1-dim Hecke eigenspace inside
//   Vtw = intersection_P [ker q_P^+(T_P) + ker q_P^-(T_P)].
// We compute, at level Pbad and Pbad*p2:
//   dim of untwisted V+ = cap ker q^+,  twisted-by-(-1) V- = cap ker q^-,  and Vtw.
// Then for each 1-dim eigenform in Vtw we list a_P and its sign eps(P) over MANY primes, count
// FAILs (a_P not a +-root), and test whether the signs form a quadratic ray-class character.
//
// Run: magma refine_serre4.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11); Rt<y> := PolynomialRing(F11);
load "mod11_trnm.txt";

LOGF := "refine_serre4_out.txt";
PrintFile(LOGF, "# refine_serre4: single-eigenform twist test" : Overwrite := true);
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

// test if signs (list of <ideal, GF2 bit>) form a quadratic char of conductor p2^a * Pbad^b * oo
function FitChar(sd, Pbad, p2, OK)
    for a in [0..6] do
      for b in [0..1] do
        m := (a eq 0 select ideal<OK|1> else p2^a) * (b eq 0 select ideal<OK|1> else Pbad);
        R, mp := RayClassGroup(m, [1,2]);
        Q, qm := quo<R | [2*R.j : j in [1..Ngens(R)]]>;
        rk := Ngens(Q);
        if rk eq 0 then
            if forall{p : p in sd | p[2] eq GF(2)!0} then return true, Sprintf("trivial (cond 1)"); end if;
            continue;
        end if;
        rows := [[GF(2)!c : c in Eltseq(qm(p[1] @@ mp))] : p in sd];
        bb := Vector(GF(2), [p[2] : p in sd]);
        A := Matrix(GF(2), rows);
        cons, sol := IsConsistent(Transpose(A), bb);
        if cons then
            nontriv := exists{c : c in Eltseq(sol) | c ne GF(2)!0};
            return true, Sprintf("conductor p2^%o * Pbad^%o, %o (func %o)", a, b,
                nontriv select "NONTRIVIAL" else "trivial", Eltseq(sol));
        end if;
      end for;
    end for;
    return false, "no 2/Pbad-power quadratic character fits";
end function;

procedure Analyze(Level, NFIND, NTEST, tag)
    LOG(Sprintf("==== level %o (norm %o) ====", tag, Norm(Level)));
    M := HilbertCuspForms(K, Level, [2,2]);
    d := Dimension(M); Id := IdentityMatrix(F11, d);
    LOG(Sprintf("  cuspidal dim = %o", d));
    Vp := VectorSpace(F11, d); Vm := Vp; Vt := Vp;
    nf := 0;
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11); T2 := T*T;
        kp := Kernel(T2 - tg[2]*T + tg[3]*Id);
        km := Kernel(T2 + tg[2]*T + tg[3]*Id);
        Vp := Vp meet kp; Vm := Vm meet km; Vt := Vt meet (kp + km);
        T := 0; T2 := 0; nf +:= 1;
        LOG(Sprintf("    [find] N=%o: dim V+=%o V-=%o Vtw=%o",
            tg[4], Dimension(Vp), Dimension(Vm), Dimension(Vt)));
        if nf ge NFIND then break; end if;
    end for;
    LOG(Sprintf("  FINAL dims: V+(untwisted)=%o  V-=%o  Vtw=%o",
        Dimension(Vp), Dimension(Vm), Dimension(Vt)));
    V := Vt;
    if Dimension(V) eq 0 then LOG("  no twisted survivor."); return; end if;
    if Dimension(V) gt 4 then LOG("  Vtw too big to read eigenforms cleanly; skipping."); return; end if;

    // try to split Vt into 1-dim eigenlines using one prime; else analyze charpoly
    BM := BasisMatrix(V); dV := Dimension(V);
    xs := [];
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        C := Solution(BM, BM*T); T := 0;
        evs := Eigenvalues(C);
        nlines := &+[Integers()| Dimension(Eigenspace(C, ev[1])) : ev in evs];
        if nlines eq dV then
            for ev in evs do for bvec in Basis(Eigenspace(C, ev[1])) do Append(~xs, bvec); end for; end for;
            LOG(Sprintf("  Vtw splits into %o eigenline(s) via N=%o", #xs, tg[4]));
            break;
        end if;
    end for;
    if #xs ne dV then
        LOG(Sprintf("  Vtw (dim %o) is NOT semisimple / not a sum of distinct eigenforms (non-diagonalizable T) -> spurious, not genuine eigenforms.", dV));
        return;
    end if;

    // per eigenline: a_P, sign, FAILs over NTEST primes
    for li in [1..#xs] do
        x := xs[li];
        okc := 0; failc := 0; ntest := 0; sd := [];
        for tg in targets do
            if Norm(tg[1] + Level) ne 1 then continue; end if;
            T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
            C := Solution(BM, BM*T); T := 0;
            w := x*C; k := 1; while x[k] eq F11!0 do k +:= 1; end while;
            a := w[k]/x[k];
            isE := (w eq a*x);
            inP := (a^2 - tg[2]*a + tg[3]) eq F11!0;
            inM := (a^2 + tg[2]*a + tg[3]) eq F11!0;
            if not isE then failc +:= 1;
            elif inP or inM then okc +:= 1;
                if tg[2] ne 0 and (inP xor inM) then
                    Append(~sd, <tg[1], inP select GF(2)!0 else GF(2)!1>);
                end if;
            else failc +:= 1; end if;
            ntest +:= 1;
            if ntest ge NTEST then break; end if;
        end for;
        LOG(Sprintf("  eigenline %o: a_P a (+-)root at %o/%o primes, FAIL %o; %o sign points",
            li, okc, ntest, failc, #sd));
        if failc eq 0 and #sd gt 0 then
            ok, desc := FitChar(sd, Pbad, p2, OK);
            LOG(Sprintf("    twist eps: %o", desc));
        end if;
    end for;
end procedure;

Analyze(Pbad,    14, 70, "Pbad");
Analyze(Pbad*p2, 14, 40, "Pbad*p2");
LOG("DONE.");
exit;
