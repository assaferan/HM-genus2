// refine_serre2.m -- CORRECT per-eigenvector refinement of the mod-11 Serre match.
//
// V (dim 2, Hecke-stable) = sum_P [ker q_P^+(T_P) + ker q_P^-(T_P)].  Its two simultaneous
// Hecke eigenlines x1,x2 are the candidate forms.  For a genuine match each eigenline's
// eigenvalue a_P must, at EVERY prime, be a root of q_P^+ (sign eps=+1) or q_P^- (eps=-1):
//     a_P^2 - tTr*a_P + tNm = 0   (+)     or     a_P^2 + tTr*a_P + tNm = 0   (-)
// and the signs eps(P) must form a single quadratic character (the twist).
//
// We split V using one prime P0 whose 2x2 restriction C_{P0} has two distinct F11
// eigenvalues; x1,x2 are then eigenvectors of every C_Q (Hecke operators commute).
//
// Run: magma refine_serre2.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
load "mod11_trnm.txt";

LOGF := "refine_serre2_out.txt";
PrintFile(LOGF, "# refine_serre2: per-eigenvector twist test" : Overwrite := true);
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

// classify eigenvalue a against q^+/q^- roots; returns "+","-","both","none"
function ClassifySign(a, tTr, tNm)
    pp := a^2 - tTr*a + tNm;
    pm := a^2 + tTr*a + tNm;
    ip := (pp eq F11!0); im := (pm eq F11!0);
    if ip and im then return "both"; elif ip then return "+"; elif im then return "-";
    else return "none"; end if;
end function;

// restrict T to subspace with basis matrix BM (BM*T = C*BM)
function Restrict(BM, T)  return Solution(BM, BM*T); end function;

// eigenvalue of 2x2 (or rxr) C on a known eigen-row x: a with x*C = a*x
function EigOn(x, C)
    w := x*C;
    k := 1; while x[k] eq F11!0 do k +:= 1; end while;
    a := w[k]/x[k];
    return a, (w eq a*x);   // scalar, is-genuine-eigenvector
end function;

procedure Analyze(Level, NMAX, tag)
    LOG(Sprintf("==== level %o (norm %o) ====", tag, Norm(Level)));
    M := HilbertCuspForms(K, Level, [2,2]);
    d := Dimension(M);
    Id := IdentityMatrix(F11, d);
    LOG(Sprintf("  cuspidal dim = %o", d));
    // PASS 1: find V
    cur := VectorSpace(F11, d); lastdim := d; stable := 0;
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        T2 := T*T;
        cur := cur meet (Kernel(T2 - tg[2]*T + tg[3]*Id) + Kernel(T2 + tg[2]*T + tg[3]*Id));
        T := 0; T2 := 0;
        if Dimension(cur) eq lastdim then stable +:= 1; else stable := 0; end if;
        lastdim := Dimension(cur);
        if Dimension(cur) eq 0 then break; end if;
        if Dimension(cur) le 3 and stable ge 4 then break; end if;
    end for;
    dV := Dimension(cur);
    LOG(Sprintf("  dim V = %o", dV));
    if dV eq 0 then LOG("  -> no survivor (clean: no match at this level)."); return; end if;
    V := cur; BM := BasisMatrix(V);

    // split V into eigenlines via a prime with distinct F11 eigenvalues on V
    xs := []; P0info := "";
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        C := Restrict(BM, T); T := 0;
        evs := Eigenvalues(C);   // set of <eval, mult>
        if #evs eq dV then       // dV distinct eigenvalues => fully splits
            for ev in evs do
                Append(~xs, Basis(Eigenspace(C, ev[1]))[1]);
            end for;
            P0info := Sprintf("split by prime N=%o into %o lines", tg[4], #xs);
            break;
        end if;
    end for;
    if #xs eq 0 then
        LOG("  V does not split into 1-dim eigenlines over F11 (irreducible/repeated) -- single eigensystem over extension.");
        return;
    end if;
    LOG("  " * P0info);

    // PASS 2: per-eigenline, classify a_P over NMAX primes
    nl := #xs;
    okcount := [0 : i in [1..nl]];
    badcount := [0 : i in [1..nl]];
    notEig := [0 : i in [1..nl]];
    signdata := [* [] : i in [1..nl] *];   // per line: list <P_ideal, GF(2) bit>
    used := 0;
    for tg in targets do
        if Norm(tg[1] + Level) ne 1 then continue; end if;
        T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
        C := Restrict(BM, T); T := 0;
        line := Sprintf("  N=%o p=%o:", tg[4], tg[5]);
        for i in [1..nl] do
            a, isE := EigOn(xs[i], C);
            if not isE then notEig[i] +:= 1; line *:= Sprintf(" [%o:notEig]", i); continue; end if;
            s := ClassifySign(a, tg[2], tg[3]);
            line *:= Sprintf(" [%o:a=%o %o]", i, a, s);
            if s eq "none" then badcount[i] +:= 1;
            else okcount[i] +:= 1;
                if tg[2] ne 0 then
                    if s eq "+" then Append(~signdata[i], <tg[1], GF(2)!0>);
                    elif s eq "-" then Append(~signdata[i], <tg[1], GF(2)!1>); end if;
                end if;
            end if;
        end for;
        LOG(line);
        used +:= 1;
        if used ge NMAX then break; end if;
    end for;
    for i in [1..nl] do
        LOG(Sprintf("  line %o: root OK at %o/%o primes, FAIL %o, notEig %o",
            i, okcount[i], used, badcount[i], notEig[i]));
    end for;

    // identify twist per line
    for i in [1..nl] do
        if badcount[i] gt 0 then
            LOG(Sprintf("  line %o: NOT a match (fails root test at %o primes).", i, badcount[i]));
            continue;
        end if;
        sd := signdata[i];
        if #sd eq 0 then LOG(Sprintf("  line %o: all tTr=0, sign undetermined.", i)); continue; end if;
        found := false;
        for e in [0..8] do
            mod_ := (e eq 0) select ideal<OK|1> else p2^e;
            R, mp := RayClassGroup(mod_, [1,2]);
            Q, qm := quo<R | [2*R.j : j in [1..Ngens(R)]]>;
            rk := Ngens(Q);
            if rk eq 0 then
                // only trivial char; consistent iff all bits 0
                if forall{p : p in sd | p[2] eq GF(2)!0} then
                    LOG(Sprintf("  line %o: eps TRIVIAL (no twist), conductor 1.", i)); found:=true; break;
                end if; continue;
            end if;
            rows := [[GF(2)!c : c in Eltseq(qm(p[1] @@ mp))] : p in sd];
            b := Vector(GF(2), [p[2] : p in sd]);
            A := Matrix(GF(2), rows);
            cons, sol := IsConsistent(Transpose(A), b);
            if cons then
                nontriv := exists{c : c in Eltseq(sol) | c ne GF(2)!0};
                LOG(Sprintf("  line %o: eps consistent at conductor p2^%o, %o (functional %o)",
                    i, e, nontriv select "NONTRIVIAL twist" else "trivial", Eltseq(sol)));
                found := true; break;
            end if;
        end for;
        if not found then LOG(Sprintf("  line %o: signs do NOT form a 2-power quadratic character.", i)); end if;
    end for;
end procedure;

Analyze(Pbad*p2, 40, "Pbad*p2");
Analyze(Pbad, 12, "Pbad");
LOG("DONE.");
exit;
