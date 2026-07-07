// kernel_sweep.m -- comprehensive Goal-1 sweep via the mod-p inert-kernel
// method (no NewformDecomposition).  Processes ALL entries; skips spaces whose
// dimension exceeds DIMCAP (dense F_p kernel infeasible) WITHOUT building Hecke
// operators (so it never OOMs).  Resume-aware; incremental output.
//
// Discrimination control: for each entry the computed Hecke matrices are reused
// to intersect kernels for a WRONG fingerprint (all c_l shifted by 1).  A real
// match => true survivor nonzero AND control survivor 0.
//
// Run: magma kernel_sweep.m
SetColumns(0);
load "goal1_data.m";
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
DIMCAP := 30000;
NINERT := 6;
INERTLMAX := 45;
RESF := "kernel_sweep_out.txt";

// resume
done := {};
if OpenTest(RESF, "r") then
    fh := Open(RESF, "r");
    while true do s := Gets(fh); if IsEof(s) then break; end if;
        if #s ge 5 and s[1..5] eq "ENTRY" then done := done join {StringToInteger(Split(s," \t")[2])}; end if;
    end while; delete fh;
else
    PrintFile(RESF, "# kernel_sweep: ENTRY idx d p realnorm -> status (trueSurv/ctrlSurv)");
end if;
procedure REC(s) printf "%o\n", s; PrintFile(RESF, s); end procedure;

function IdOfNorm(OF, N)
    J := ideal<OF|1>;
    for pe in Factorization(N) do
        P := Factorization(pe[1]*OF)[1][1];
        J := J * P^(pe[2] div Valuation(Norm(P), pe[1]));
    end for;
    return J;
end function;

for ent in entries do
    idx := ent[1]; d := ent[2]; p := ent[3]; serre := ent[4]; fc := ent[5]; hc := ent[6]; cond := ent[7];
    if idx in done then continue; end if;
    F<w> := QuadraticField(d); OF := Integers(F); Fp := GF(p);
    C := HyperellipticCurve(&+[fc[i]*x^(i-1):i in [1..#fc]], &+[hc[i]*x^(i-1):i in [1..#hc]]);
    dC := Integers()!Discriminant(C);
    lev := IdOfNorm(OF, serre);
    if cond mod p eq 0 then lev := lev * Factorization(p*OF)[1][1]^2; end if;
    hdr := Sprintf("ENTRY %o d%o p%o realnorm%o", idx, d, p, Norm(lev));

    M := HilbertCuspForms(F, lev, [2,2]);
    dm := Dimension(M);
    if dm eq 0 then REC(Sprintf("%o -> EMPTY", hdr)); continue; end if;
    if dm gt DIMCAP then REC(Sprintf("%o -> SKIP-BIG dim %o", hdr, dm)); continue; end if;

    // inert conditions from A: c_l = -a2(A,l) mod p
    conds := [];
    for l in PrimesInInterval(3, 200) do
        if #conds ge NINERT then break; end if;
        if l eq p or (d mod l eq 0) or (dC mod l eq 0) or (l gt INERTLMAX) then continue; end if;
        if IsSquare(GF(l)!d) then continue; end if;
        if Norm(lev) mod l eq 0 then continue; end if;
        c := Fp!(-(Integers()!Coefficient(LPolynomial(ChangeRing(C,GF(l))), 2)));
        Append(~conds, <Factorization(l*OF)[1][1], c, l>);
    end for;

    t0 := Cputime();
    V := VectorSpace(Fp, dm); I := IdentityMatrix(Fp, dm);
    Ts := []; cs := []; ls := [];
    for cc in conds do
        Append(~Ts, ChangeRing(Matrix(HeckeOperator(M, cc[1])), Fp));
        Append(~cs, cc[2]); Append(~ls, cc[3]);
    end for;
    // true fingerprint
    cur := V; for j in [1..#Ts] do cur := cur meet Kernel(Ts[j] - cs[j]*I); if Dimension(cur) eq 0 then break; end if; end for;
    trueS := Dimension(cur);
    // control: shift every c by 1
    cur2 := V; for j in [1..#Ts] do cur2 := cur2 meet Kernel(Ts[j] - (cs[j]+1)*I); if Dimension(cur2) eq 0 then break; end if; end for;
    ctrlS := Dimension(cur2);

    status := (trueS gt 0 and ctrlS eq 0) select "MATCH" else
              ((trueS gt 0 and ctrlS gt 0) select "MATCH?(ctrl nonzero)" else "NO-MATCH");
    REC(Sprintf("%o -> %o  (true %o / ctrl %o)  dim %o  inert %o  [%.1os]",
        hdr, status, trueS, ctrlS, dm, ls, Cputime(t0)));
end for;
REC("# kernel_sweep pass done");
exit;
