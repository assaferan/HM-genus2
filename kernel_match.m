// kernel_match.m -- LIGHTER method: find the matching Hilbert eigenform mod p
// WITHOUT NewformDecomposition.  On the full cusp space M = HilbertCuspForms,
// intersect ker(T_P - c_l I) mod p over INERT primes l, where the induced
// structure forces the Hecke eigenvalue a_P(f) = -a2(A,l) mod p  (a2 = T^2
// coeff of A's L-polynomial at l).  At inert primes this value lies in F_p even
// though split-prime eigenvalues lie in F_{p^2}, so the kernels are over F_p.
//
// Skips char-0 decomposition (the memory/time bottleneck), so it reaches larger
// levels.  Reports dim of the surviving mod-p Hecke eigenspace (>0 = match).
//
// Usage: set IDXS below.  Run: magma kernel_match.m

SetColumns(0);
load "goal1_data.m";
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
IDXS := [12, 13, 16, 17];   // mod-5 mid-range (serre 144/162 -> norm 90000/101250)
NINERT := 10;               // inert primes to intersect
INERTLMAX := 40;            // cap l so Hecke prime norm l^2 stays cheap
DIMCAP := 40000;            // skip Hecke/kernel if full-space dim exceeds this

function IdOfNorm(OF, N)     // one ideal of prime-to-p norm N (squarefree-ish)
    J := ideal<OF|1>;
    for pe in Factorization(N) do
        P := Factorization(pe[1]*OF)[1][1];
        fdeg := Valuation(Norm(P), pe[1]);
        J := J * P^(pe[2] div fdeg);
    end for;
    return J;
end function;

LOGF := "kernel_match_out.txt";
PrintFile(LOGF, "# kernel_match: mod-p inert kernel intersection" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

for target in IDXS do
    ent := [e : e in entries | e[1] eq target][1];
    idx := ent[1]; d := ent[2]; p := ent[3]; serre := ent[4]; fc := ent[5]; hc := ent[6]; cond := ent[7];
    F<w> := QuadraticField(d); OF := Integers(F); Fp := GF(p);
    C := HyperellipticCurve(&+[fc[i]*x^(i-1):i in [1..#fc]], &+[hc[i]*x^(i-1):i in [1..#hc]]);
    dC := Integers()!Discriminant(C);

    // build level: SerreIdeal (* Pp^2 if p is bad)
    lev := IdOfNorm(OF, serre);
    pbad := (cond mod p eq 0);
    if pbad then lev := lev * Factorization(p*OF)[1][1]^2; end if;
    LOG(Sprintf("\n===== idx %o  d%o p%o  level norm %o (pbad=%o) =====", idx, d, p, Norm(lev), pbad));

    // inert-prime conditions from A: <ideal (l), c_l = -a2(A,l) mod p>
    conds := [];
    for l in PrimesInInterval(3, 200) do
        if #conds ge NINERT then break; end if;
        if l eq p or (d mod l eq 0) or (dC mod l eq 0) or (l gt INERTLMAX) then continue; end if;
        if IsSquare(GF(l)!d) then continue; end if;          // want INERT
        if Norm(lev) mod l eq 0 then continue; end if;
        LA := LPolynomial(ChangeRing(C, GF(l)));
        c := Fp ! (-(Integers()!Coefficient(LA, 2)));
        Pin := Factorization(l*OF)[1][1];    // = (l), norm l^2
        Append(~conds, <Pin, c, l>);
    end for;
    LOG(Sprintf("  %o inert conditions (l): %o", #conds, [x[3]: x in conds]));

    t0 := Cputime();
    M := HilbertCuspForms(F, lev, [2,2]);
    dm := Dimension(M);
    LOG(Sprintf("  built space: dim %o  (%.1o s)", dm, Cputime(t0)));
    if dm eq 0 then LOG("  empty space"); continue; end if;
    if dm gt DIMCAP then
        LOG(Sprintf("  dim %o > DIMCAP %o -- dense F_%o kernel infeasible, SKIP", dm, DIMCAP, p));
        continue;
    end if;

    V := VectorSpace(Fp, dm); cur := V; I := IdentityMatrix(Fp, dm);
    for cc in conds do
        t1 := Cputime();
        Tp := ChangeRing(Matrix(HeckeOperator(M, cc[1])), Fp);
        ker := Kernel(Tp - cc[2]*I);
        cur := cur meet ker;
        LOG(Sprintf("    l=%o (c=%o): kerdim %o  cumulative %o   [Hecke %.1os]",
            cc[3], cc[2], Dimension(ker), Dimension(cur), Cputime(t1)));
        if Dimension(cur) eq 0 then break; end if;
    end for;
    LOG(Sprintf("  >>> surviving mod-%o eigenspace dim = %o   (%o => MATCH)", p, Dimension(cur),
        Dimension(cur) gt 0 select "nonzero" else "zero"));
    LOG(Sprintf("  total time %.1os", Cputime(t0)));
end for;
LOG("\nDONE.");
exit;
