// cfe_conductor2.m -- determine c = conductor exponent of A_true at 2 via the functional equation.
// Precompute A_true's Euler factors up to the cutoff, then search c in N = 5^4*66179*4^c.
//
// Run: magma cfe_conductor2.m   (background; ~minutes for the precompute)

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K);
ZT<T> := PolynomialRing(Integers());

LOGF := "cfe_conductor2_out.txt";
PrintFile(LOGF, "# cfe_conductor2" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K); Cint := IntegralModelOf(Cmin);
f, h := HyperellipticPolynomials(Cint);
F0 := PK![OK!c : c in Coefficients(f + (Parent(f)!h)^2/4)];
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(&+[ideal<OK|c> : c in Coefficients(F0)|c ne 0])]);
_, alpha := IsPrincipal(Msq);
F1 := PK![ (OK!c)/alpha^2 : c in Coefficients(F0) ];
_, g := IsPrincipal(&+[ideal<OK|c> : c in Coefficients(F1)|c ne 0]);
Fprim := PK![ c/g : c in Coefficients(F1) ];
LOG(Sprintf("setup done; g norm %o", Norm(g)));

function ChiG(P)
    Fp, red := ResidueClassField(P); v := red(OK!g);
    if v eq 0 then return 0; end if; return IsSquare(v) select 1 else -1;
end function;

// reduce a (possibly non-integral) sextic model at P; return LPolynomial or 0 if not found
function ReduceL(Fpoly, P)
    Fp, redP := ResidueClassField(P); Rx<xx> := PolynomialRing(Fp); n := 6;
    vals := [Valuation(c, P) : c in Coefficients(Fpoly) | c ne 0];
    m := Min(vals); pi := UniformizingElement(P);
    Fs := PK![ c/pi^m : c in Coefficients(Fpoly) ];
    shifts := [OK| sh : sh in [0..12] ] cat [OK| sh + sh2*s5 : sh in [0..3], sh2 in [1..3] ];
    for flip in [0,1] do
        base := flip eq 0 select Fs else PK![Coefficient(Fs, n-i) : i in [0..n]];
        for sh in shifts do
            fp := Rx![redP(OK!c) : c in Coefficients(PK!Evaluate(base, x + sh))];
            if Degree(fp) eq n and IsSquarefree(fp) then return LPolynomial(HyperellipticCurve(fp)); end if;
        end for;
    end for;
    return ZT!0;
end function;

// frobdata (A_recon good, norm<=2000): <N,p,rid,[c0..c4]>
load "mod11_frobdata.txt";
frob := AssociativeArray();
for row in frobdata do frob[<row[2], row[3]>] := &+[ZT| row[4][i+1]*T^i : i in [0..4]]; end for;

CUTOFF := 47000;
reconbad := {499, 13711};   // recon-bad primes below cutoff (A_true good); 88301,231481,66179 > cutoff
Ltrue := AssociativeArray();
nfail := 0;
for p in PrimesUpTo(CUTOFF) do
    if p eq 2 then continue; end if;
    for pe in Factorization(p*OK) do
        P := pe[1]; fP := Degree(ResidueClassField(P));
        if fP eq 1 then _rf,_rd := ResidueClassField(P); rid := Integers()!_rd(OK!s5); else rid := -1; end if;
        cg := ChiG(P);
        Lrec := ZT!0;
        if p in reconbad then
            Lrec := ReduceL(Fprim, P);          // A_true directly (good here) -- already "true", no twist
            if Lrec ne 0 then Ltrue[<p,rid>] := Lrec; continue; end if;
        end if;
        if IsDefined(frob, <p,rid>) then
            Lrec := frob[<p,rid>];
        else
            Lrec := ReduceL(F1, P);             // A_recon
        end if;
        if Lrec eq 0 then nfail +:= 1; LOG(Sprintf("  FAIL reduce p=%o rid=%o norm=%o", p, rid, Norm(P))); continue; end if;
        Ltrue[<p,rid>] := Evaluate(Lrec, cg*T);   // twist to A_true
    end for;
end for;
LOG(Sprintf("precompute done: %o factors, %o failures", #Keys(Ltrue), nfail));

function LocalFactor(p, deg)
    if p eq 2 then return ZT!1; end if;
    res := ZT!1;
    for pe in Factorization(p*OK) do
        P := pe[1]; fP := Degree(ResidueClassField(P));
        if fP eq 1 then _rf,_rd := ResidueClassField(P); rid := Integers()!_rd(OK!s5); else rid := -1; end if;
        if not IsDefined(Ltrue, <p,rid>) then return ZT!1; end if;   // beyond precompute (shouldn't happen < cutoff)
        res := res * Evaluate(Ltrue[<p,rid>], T^fP);
    end for;
    return res mod T^(deg+1);
end function;
cf := func< p, d | LocalFactor(p, d) >;

for c in [6..11] do
    N := 5^4 * 66179 * 4^c;
    L := LSeries(2, [0,0,0,0,1,1,1,1], N, cf : Sign := 0, Precision := 1);
    err := 999.0;
    try err := CheckFunctionalEquation(L); catch e; LOG(Sprintf("  c=%o: CFE err %o", c, e`Object)); continue; end try;
    LOG(Sprintf("  c=%o (N=%o): |CFE| = %o", c, N, err));
end for;
LOG("DONE.");
exit;
