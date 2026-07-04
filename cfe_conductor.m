// cfe_conductor.m -- determine the conductor exponent c = a_{p2}(A_true) at the inert prime 2
// via the functional equation of L(A_true/K, s) = L(W/Q, s), W = Res_{K/Q} A_true.
// A_true = A_recon (x) chi_g is good at 499,13711,88301,231481,huge; conductor = p * (2)^c,
// so N(W) = disc_K^4 * N_{K/Q}(cond) = 5^4 * 66179 * 4^c.  Search c via CheckFunctionalEquation.
//
// Euler factors: A_true's L_P = chi_g-twist of A_recon's (from frobdata): L_P^true(T)=L_P^rec(chi_g(P)T).
// Degree 8 over Q, motivic weight 1 => gamma = [0,0,0,0,1,1,1,1], weight 2. Factor at 2 = 1 (additive).
//
// Run: magma cfe_conductor.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K);
ZT<T> := PolynomialRing(Integers());

LOGF := "cfe_conductor_out.txt";
PrintFile(LOGF, "# cfe_conductor" : Overwrite := true);
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
LOG(Sprintf("g norm %o", Norm(g)));

// frobdata: A_recon good-prime L-polys  <N,p,rid,[c0..c4]>
load "mod11_frobdata.txt";
frob := AssociativeArray();
for row in frobdata do frob[<row[2], row[3]>] := row[4]; end for;

function ChiG(P)
    Fp, red := ResidueClassField(P); v := red(OK!g);
    if v eq 0 then return 0; end if;
    return IsSquare(v) select 1 else -1;
end function;

// local factor at rational prime p, as a polynomial in T (=p^{-s}), truncated to degree deg
function LocalFactor(p, deg)
    if p eq 2 then return ZT!1; end if;          // additive at 2
    res := ZT!1;
    for pe in Factorization(p*OK) do
        P := pe[1]; fP := Degree(ResidueClassField(P));   // residue degree (1 split, 2 inert)
        // recover rid = red(s5) for split primes to match frobdata key
        if fP eq 1 then
            Fp, red := ResidueClassField(P);
            rid := Integers()!red(OK!s5);
        else
            rid := -1;
        end if;
        key := <p, rid>;
        if not IsDefined(frob, key) then
            LOG(Sprintf("    MISSING factor: p=%o rid=%o norm=%o", p, rid, Norm(P)));
            return ZT!1;   // placeholder to avoid crash while diagnosing the cutoff
        end if;
        cs := frob[key];                       // [c0..c4], L_rec(T)=sum c_i T^i
        cg := ChiG(P);
        Lrec := &+[ ZT| cs[i+1]*T^i : i in [0..4] ];
        Ltrue := Evaluate(Lrec, cg*T);         // twist
        // substitute T -> T^fP (Norm(P)=p^fP), i.e. local factor in x=p^{-s}
        loc := Evaluate(Ltrue, T^fP);
        res := res*loc;
    end for;
    return res mod T^(deg+1);
end function;

cf := func< p, d | LocalFactor(p, d) >;

for c in [6] do
    N := 5^4 * 66179 * 4^c;
    L := LSeries(2, [0,0,0,0,1,1,1,1], N, cf : Sign := 0, Precision := 1);
    err := -1.0;
    try err := CheckFunctionalEquation(L); catch e; LOG(Sprintf("  c=%o: CFE error %o", c, e`Object)); continue; end try;
    LOG(Sprintf("  c=%o (N=%o): |CFE| = %o", c, N, err));
end for;
LOG("DONE.");
exit;
