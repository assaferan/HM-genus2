// local_at_11.m -- reduction type of Curve A at the two primes above 11 (11 splits in Q(sqrt5)).
// Determines ordinary vs supersingular (p-rank) of A = Jac(C) at each l|11, which constrains
// the Serre weight of the 2-dim pieces rho_i.
//
// Run: magma local_at_11.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, PrimesAbove;
K<s5> := QuadraticField(5); OK := Integers(K);

LOGF := "local_at_11_out.txt";
PrintFile(LOGF, "# local_at_11: reduction type of A at primes above 11" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787),
    1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089),
    1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K);
Cmin, dtw := MinimalTwist(C, K);
LOG(Sprintf("minimal twist factor %o", dtw));

function MakeIntegral(CC)
    f0, h0 := HyperellipticPolynomials(CC);
    m := LCM([Integers()| Denominator(c) : c in Coefficients(f0) cat Coefficients(h0) cat [BaseRing(CC)!1]]);
    return HyperellipticCurve(Parent(f0)!(m^2*f0), Parent(f0)!(m*h0));
end function;
Cint := MakeIntegral(Cmin);
f0, h0 := HyperellipticPolynomials(Cint);

primes11 := PrimesAbove(K, 11);
LOG(Sprintf("primes above 11: norms %o", [Norm(P) : P in primes11]));

for P in primes11 do
    Fq, red := ResidueClassField(P);
    LOG(Sprintf("== prime of norm %o, residue field size %o ==", Norm(P), #Fq));
    Rq := PolynomialRing(Fq);
    fb := Rq![red(c) : c in Coefficients(f0)];
    hb := Rq![red(c) : c in Coefficients(h0)];
    try
        Cbar := HyperellipticCurve(fb, hb);
        if not IsNonsingular(Cbar) then LOG("  SINGULAR reduction (bad model here)"); continue; end if;
        L := LPolynomial(Cbar);     // T^4 - a1 T^3 + a2 T^2 - p a1 T + p^2
        LOG(Sprintf("  L-poly: %o", L));
        cs := Coefficients(L);      // [c0,c1,c2,c3,c4] = [p^2, -p a1, a2, -a1, 1]
        a1 := -cs[4]; a2 := cs[3];
        LOG(Sprintf("  a1 = %o (mod 11 = %o), a2 = %o (mod 11 = %o)", a1, a1 mod 11, a2, a2 mod 11));
        // Newton polygon / p-rank at 11
        NP := NewtonPolygon(L, 11);
        slopes := Slopes(NP);
        prank := #[s : s in slopes | s eq 0];
        LOG(Sprintf("  Newton slopes at 11: %o  => p-rank = %o", slopes, prank));
        if prank eq 2 then LOG("  => ORDINARY at this prime (p-rank 2)");
        elif prank eq 0 then LOG("  => SUPERSINGULAR (p-rank 0)");
        else LOG(Sprintf("  => p-rank %o (mixed)", prank)); end if;
        // also factor L mod 11 to see the local mod-11 structure
        L11 := ChangeRing(L, GF(11));
        LOG(Sprintf("  L mod 11 factors: %o", Factorization(L11)));
    catch e
        LOG(Sprintf("  error reducing: %o", e`Object));
    end try;
end for;
LOG("DONE.");
exit;
