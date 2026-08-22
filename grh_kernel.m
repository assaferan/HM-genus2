// grh_kernel.m -- Step 2: GRH-conditional trace-agreement certificate sigma = rho_f for a
// torsion-type curve, sourcing a_P(f) mod lambda from the kernel SURVIVOR EIGENVECTOR v
// (validated in grh_kernel_validate.m) instead of a char-0 NewformDecomposition. This reaches
// levels far beyond decomposition.
//
//   tr sigma(Frob_P) = -#C(F_P)            (mod l)          [from the curve]
//   tr rho_f(Frob_P) = c_P, where v*T_P = c_P*v over F_l    [from the survivor v]
//
// sigma, rho_f are 2-dim, det = chi_l (torsion structure), unramified outside {level, l}; we
// check sigma irreducible (a Frobenius char poly irreducible over F_l) and verify EXACT trace
// agreement for all good primes N(P) <= BOUND (>> the GRH Faltings-Serre bound ~ (log cond)^2).
//
// Usage:  magma idx:=3 grh_kernel.m            (idx into torsion_data; optional BOUND:=... )
// Writes grh_kernel_<idx>.out. Sentinel on success: "GRH-KERNEL CERT <label>: MODULAR".
SetColumns(0);
load "torsion_data.m";
if not assigned BOUND then BOUND := 3000; end if;
BOUND := StringToInteger(Sprintf("%o", BOUND));
idx := StringToInteger(idx);
row := torsion_data[idx];
d := row[1]; l := row[2]; lab := row[3]; fc := row[4]; hc := row[5]; N := row[6]; cc := row[7];

outf := Sprintf("grh_kernel_%o.out", idx);
PrintFile(outf, Sprintf("# grh_kernel idx %o: %o (d=%o l=%o), BOUND=%o", idx, lab, d, l, BOUND) : Overwrite := true);
procedure LOG(outf, s) printf "%o\n", s; PrintFile(outf, s); end procedure;

K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K); Fl := GF(l); RT := PolynomialRing(Fl);
f := R![K| c[1]+c[2]*a : c in fc]; hh := R![K| c[1]+c[2]*a : c in hc];
C := HyperellipticCurve(f, hh);
p2 := Factorization(2*OK)[1][1];
if d eq 2 then base := ideal<OK | cc[1]+cc[2]*a>;
else base := 1*OK; for pe in Factorization(Conductor(C)) do if Norm(pe[1]) mod 2 ne 0 then base := base*pe[1]^pe[2]; end if; end for;
end if;
emax := Valuation(Conductor(C), p2);

function tP(P)
    k,red := ResidueClassField(P); Rk := PolynomialRing(k);
    Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(f)], Rk![red(c):c in Coefficients(hh)]);
    return Fl!(-#Ck);
end function;

// fingerprint primes for isolating the survivor
fpP := []; fpT := [];
for pp in PrimesUpTo(120) do
    if #fpP ge 14 then break; end if;
    if pp eq 2 or pp eq l then continue; end if;
    for tup in Factorization(pp*OK) do
        P := tup[1]; if Norm(base+P) ne 1 then continue; end if;
        Append(~fpP, P); Append(~fpT, tP(P));
        if #fpP ge 14 then break; end if;
    end for;
end for;

// ---- find the match level (smallest e: survivor 1-dim, control collapses) and get v ----
Mmatch := 0; vv := 0; jj := 0; ematch := -1; Nlev := 0;
for e in [0..emax] do
    Nl := base * p2^e;
    M := HilbertCuspForms(K, Nl, [2,2]); dm := Dimension(M);
    if dm eq 0 then continue; end if;
    V := VectorSpace(Fl, dm); I := IdentityMatrix(Fl, dm);
    // survivor under true fingerprint
    cur := V; for i in [1..#fpP] do
        Tp := ChangeRing(Matrix(HeckeOperator(M, fpP[i])), Fl);
        cur := cur meet Kernel(Tp - fpT[i]*I); if Dimension(cur) eq 0 then break; end if;
    end for;
    if Dimension(cur) ne 1 then continue; end if;
    // control t->t+1 must collapse
    cc2 := V; for i in [1..#fpP] do
        Tp := ChangeRing(Matrix(HeckeOperator(M, fpP[i])), Fl);
        cc2 := cc2 meet Kernel(Tp - (fpT[i]+1)*I); if Dimension(cc2) eq 0 then break; end if;
    end for;
    if Dimension(cc2) ne 0 then continue; end if;
    Mmatch := M; vv := Basis(cur)[1]; jj := rep{ i : i in [1..dm] | Basis(cur)[1][i] ne 0 };
    ematch := e; Nlev := Nl;
    LOG(outf, Sprintf("  match level: e=%o levelN=%o dim=%o survivor=1 control=0", e, Norm(Nl), dm));
    break;
end for;
error if ematch lt 0, "no match level found";

// eigenvalue of T_Q on the survivor:  v * T_Q = c * v
function survEig(Q)
    Tq := ChangeRing(Matrix(HeckeOperator(Mmatch, Q)), Fl);
    w := vv*Tq; c := w[jj]/vv[jj];
    error if w ne c*vv, Sprintf("survivor not an eigenvector at N(Q)=%o", Norm(Q));
    return c;
end function;

// ---- sigma irreducible? (a Frobenius char poly irreducible over F_l) ----
irr := false;
for i in [1..#fpP] do
    cp := RT ! (RT.1^2 - (Integers()!fpT[i])*RT.1 + (Norm(fpP[i]) mod l));
    if IsIrreducible(cp) then irr := true; break; end if;
end for;
LOG(outf, Sprintf("  sigma irreducible: %o", irr));

// ---- verify tr sigma = tr rho_f for all good P up to BOUND ----
LOG(outf, Sprintf("  verifying trace agreement, good primes N(P) <= %o ...", BOUND));
// Progress step must scale with BOUND: a fixed 1000 meant that every run with BOUND <= 1000
// (i.e. every large-dim shard, which is exactly where progress matters) printed NOTHING between
// "verifying ..." and the final total -- 12+ hours of silence with no way to tell how far along.
logstep := Max(50, BOUND div 10);
nchk := 0; fails := []; t0 := Cputime(); nextlog := logstep;
for pp in PrimesUpTo(BOUND) do
    if pp eq 2 or pp eq l then continue; end if;
    for tup in Factorization(pp*OK) do
        P := tup[1];
        if Norm(P) gt BOUND or Norm(Nlev+P) ne 1 then continue; end if;
        if tP(P) ne survEig(P) then Append(~fails, Norm(P)); end if;
        nchk +:= 1;
    end for;
    if pp gt nextlog then LOG(outf, Sprintf("    ... up to %o: %o primes, %o disagreements (%.1o s)", pp, nchk, #fails, Cputime(t0))); nextlog +:= logstep; end if;
end for;

LOG(outf, Sprintf("  TOTAL: %o good primes N(P) <= %o; DISAGREEMENTS: %o", nchk, BOUND, #fails));
if #fails eq 0 and irr then
    LOG(outf, Sprintf("%-11o d=%o l=%o levelN=%o e=%o primes=%o disagree=0 irred  GRH-KERNEL CERT %o: MODULAR",
        lab, d, l, Norm(Nlev), ematch, nchk, lab));
else
    LOG(outf, Sprintf("%-11o  GRH-KERNEL CERT %o: INCOMPLETE (disagree=%o irr=%o)", lab, lab, #fails, irr));
end if;
LOG(outf, Sprintf("GRH-KERNEL DONE %o", idx));
exit;
