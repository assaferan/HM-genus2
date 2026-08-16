// kernel_torsion.m -- match sigma (2-dim sub of A[l], torsion-type) to a Hilbert
// newform WITHOUT NewformDecomposition, via mod-l Hecke eigen-kernel intersection.
//
// On the full cusp space M = HilbertCuspForms(F, base*p2^e, [2,2]) intersect
//     survivor = MEET_P ker(T_P - t_P I)   over F_l,
// where t_P = -#C(F_P) mod l is the torsion fingerprint (all good primes P, split
// AND inert -- sigma is a genuine GL_2/F rep, so t_P in F_l directly).  survivor dim
// > 0  <=>  a newform matches.  Discrimination control: shift t_P -> t_P + 1 must
// collapse survivor to 0.  Sweeps the 2-part exponent e; first e with (true>0 and
// control=0) is the match (mod-l level-lowering minimal level).
//
// Usage:  magma idx:=3 kernel_torsion.m        (one curve, 1-based index into torsion_data)
// Optional: NCOND:=... (default 14 good primes), PBOUND:=... (default 120).

SetColumns(0);
load "torsion_data.m";
if not assigned NCOND then NCOND := 14; end if;
if not assigned PBOUND then PBOUND := 120; end if;
NCOND := StringToInteger(Sprintf("%o", NCOND));
PBOUND := StringToInteger(Sprintf("%o", PBOUND));
idx := StringToInteger(idx);
row := torsion_data[idx];
d := row[1]; l := row[2]; lab := row[3]; fc := row[4]; hc := row[5]; N := row[6]; cc := row[7];

outf := Sprintf("kernel_%o.out", idx);
PrintFile(outf, Sprintf("# kernel_torsion shard %o: %o (d=%o l=%o)", idx, lab, d, l) : Overwrite := true);
procedure LOG(outf, s) printf "%o\n", s; PrintFile(outf, s); end procedure;

K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K); Fl := GF(l);
f := R![K| c[1]+c[2]*a : c in fc]; hh := R![K| c[1]+c[2]*a : c in hc];
C := HyperellipticCurve(f, hh);
p2 := Factorization(2*OK)[1][1];
if d eq 2 then base := ideal<OK | cc[1]+cc[2]*a>;
else
    Cnd := Conductor(C); base := 1*OK;
    for pe in Factorization(Cnd) do if Norm(pe[1]) mod 2 ne 0 then base := base*pe[1]^pe[2]; end if; end for;
end if;
emax := Valuation(Conductor(C), p2);
LOG(outf, Sprintf("  base-cond norm %o, e in [0..%o]", Norm(base), emax));

// --- fingerprint: t_P = -#C(F_P) mod l at good primes P coprime to base, 2, l ---
fp := [];
for pp in PrimesUpTo(PBOUND) do
    if #fp ge NCOND then break; end if;
    if pp eq 2 or pp eq l then continue; end if;
    for tup in Factorization(pp*OK) do
        P := tup[1]; if Norm(base+P) ne 1 then continue; end if;
        k,red := ResidueClassField(P); Rk := PolynomialRing(k); ok := true; Ck := 0;
        try Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(f)], Rk![red(c):c in Coefficients(hh)]);
        catch e; ok := false; end try;
        if ok then Append(~fp, <P, Fl!(-#Ck)>); end if;
        if #fp ge NCOND then break; end if;
    end for;
end for;
LOG(outf, Sprintf("  %o fingerprint primes (norms): %o", #fp, [Norm(pr[1]) : pr in fp]));

// --- intersect kernels over F_l; return surviving eigenspace dim (tw shifts t_P) ---
function Survivor(M, dm, fp, Fl, tw)
    V := VectorSpace(Fl, dm); cur := V; I := IdentityMatrix(Fl, dm);
    for pr in fp do
        Tp := ChangeRing(Matrix(HeckeOperator(M, pr[1])), Fl);
        cur := cur meet Kernel(Tp - (pr[2]+tw)*I);
        if Dimension(cur) eq 0 then break; end if;
    end for;
    return Dimension(cur);
end function;

matched := false;
for e in [0..emax] do
    Nlev := base * p2^e;
    t0 := Cputime();
    M := HilbertCuspForms(K, Nlev, [2,2]); dm := Dimension(M);
    if dm eq 0 then LOG(outf, Sprintf("  e=%o levelN=%o dim=0 (skip)", e, Norm(Nlev))); continue; end if;
    surv := Survivor(M, dm, fp, Fl, 0);
    ctrl := surv gt 0 select Survivor(M, dm, fp, Fl, 1) else 0;
    tag := (surv gt 0 and ctrl eq 0) select "MATCH" else (surv gt 0 select "CONTROL-FAIL" else "no-survivor");
    LOG(outf, Sprintf("  e=%o levelN=%o dim=%o  survivor=%o control=%o  %o  [%.1os]",
        e, Norm(Nlev), dm, surv, ctrl, tag, Cputime(t0)));
    if surv gt 0 and ctrl eq 0 then
        LOG(outf, Sprintf("%-11o d=%o l=%o levelN=%o dim=%o survivor=%o control=0  MATCH (e=%o)", lab,d,l,Norm(Nlev),dm,surv,e));
        matched := true; break;
    end if;
end for;
if not matched then LOG(outf, Sprintf("%-11o d=%o l=%o  NO MATCH up to e=%o", lab,d,l,emax)); end if;
LOG(outf, Sprintf("KERNEL SHARD %o DONE", idx));
exit;
