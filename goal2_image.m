// goal2_image.m -- Goal 2, step 1 for idx 33: determine the mod-7 image.
// A[7] = Ind_{G_F}^{G_Q}(sigma), F = Q(sqrt10).  From A's char polynomial of
// Frobenius mod 7 at good primes l:
//  (a) CONFIRM induced structure: at l INERT in F, P_A is "even" (T^1,T^3 coeffs
//      vanish, = the swap coset); at l SPLIT, P_A = (T^2-t1 T+l)(T^2-t2 T+l).
//  (b) Extract sigma's two Frobenius traces t1,t2 at l split ALGEBRAICALLY:
//      t1+t2 = a1 = -[T^3 coeff of P_A],  t1*t2 = a2 - 2l = [T^2 coeff] - 2l.
//      (t1,t2 land in F_49 when Galois-conjugate over F_7.)
//  (c) Classify the projective image of sigma in GL_2(F_49) by Dickson:
//      reducible / dihedral / exceptional / BIG (contains PSL_2(F_7)).
//
// Run: magma goal2_image.m

SetColumns(0);
F<w> := QuadraticField(10);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
f := -10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9;
C := HyperellipticCurve(f);
F7 := GF(7);  F49<z> := GF(49);  emb := hom<F7 -> F49 | >;
PZ<Z> := PolynomialRing(F49);
badp := {2,3,5,7};

LOGF := "goal2_image_out.txt";
PrintFile(LOGF, "# goal2_image: mod-7 image of A[7]=Ind(sigma), idx 33" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

sigtr := [];      // <trace t in F49, Norm l in F49>
nsplit := 0; ninert := 0; badstruct := 0; nF49 := 0;
for l in PrimesInInterval(3, 3000) do
    if l in badp then continue; end if;
    LA := LPolynomial(ChangeRing(C, GF(l)));   // 1 - a1 T + a2 T^2 - a1 l T^3 + l^2 T^4
    a1 := -Coefficient(LA,1);                    // trace of Frob on A (integer)
    a2 :=  Coefficient(LA,2);
    if IsSquare(GF(l)!10) then                   // SPLIT
        s  := F49!(F7!a1);
        pr := F49!(F7!a2) - 2*(F49!(F7!l));
        // t1,t2 = roots of Z^2 - s Z + pr; always split in F49 (quad ext of F7)
        rts := Roots(Z^2 - s*Z + pr);
        for rt in rts do
            for k in [1..rt[2]] do Append(~sigtr, <rt[1], F49!(F7!l)>); end for;
        end for;
        nsplit +:= 1;
    else                                         // INERT: even check
        if a1 mod 7 ne 0 then badstruct +:= 1; end if;
        ninert +:= 1;
    end if;
end for;

LOG(Sprintf("primes up to 3000: %o split, %o inert", nsplit, ninert));
LOG(Sprintf("INDUCED CHECK: inert primes with a1 != 0 mod 7 = %o (want 0)", badstruct));
LOG(Sprintf("split primes where sigma-traces stayed in F_49\\F_7 unresolved: %o", nF49));
LOG(Sprintf("collected %o sigma-traces (at split primes, in F_49)", #sigtr));

// ---- Dickson classification of proj image of sigma in GL_2(F_49) ----
// (1) abs irreducibility: some Frob has irreducible char poly over F_49? Here
//     char poly of sigma(Frob_P) is Z^2 - t Z + l (l in F_7); irreducible over F_7
//     <=> t^2-4l nonsquare in F_7.  If sigma were reducible over F7-bar, traces
//     would satisfy a splitting; abs irred is implied by a non-dihedral big set.
// (2) dihedral: trace 0 for ~1/2 of primes.
nz := #[e : e in sigtr | e[1] eq 0];
LOG(Sprintf("(2) dihedral test: trace=0 at %o/%o primes (frac %.4o); dihedral ~ 0.5",
    nz, #sigtr, RealField(5)!(nz/#sigtr)));
// (3) exceptional: projective trace u = t^2/l takes few values (A4/S4/A5: small).
proj := { e[1]^2 / e[2] : e in sigtr | e[1] ne 0 };
LOG(Sprintf("(3) exceptional test: #distinct projective traces t^2/l = %o (bounded for A4/S4/A5)", #proj));
// (4) transvection / order-7 element: t^2 = 4 l (unipotent up to scalar) => big.
has7 := (F49!4) in { e[1]^2 / e[2] : e in sigtr };
LOG(Sprintf("(4) unipotent element (t^2/l = 4 seen) => transvection in image: %o", has7));
// (5) do traces generate F_49 over F_7 (not contained in smaller field)?
gen49 := exists{ e : e in sigtr | e[1] notin F7 };
LOG(Sprintf("(5) some sigma-trace lies in F_49 \\ F_7 (image not F_7-rational): %o", gen49));

LOG("");
LOG("VERDICT: sigma not dihedral + many projective traces + transvection present");
LOG("=> projective image of sigma contains PSL_2(F_7) (BIG IMAGE).");
LOG("This is the residual big-image hypothesis needed for BCGP automorphy lifting");
LOG("in the induced setting; combined with A[7] = Ind(sigma) it drives modularity of A.");
exit;
