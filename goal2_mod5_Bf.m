// goal2_mod5_Bf.m -- gather the data of the newform f (idx 5, orbit 6) needed to
// pin down / construct its attached abelian surface B_f over F = Q(sqrt2).
//
// f has weight [2,2], level 𝔫 = P2^3*P5^2 (norm 5000), degree-2 Hecke field K_f.
// So B_f is an abelian surface over F of GL_2-type with RM by (an order in) K_f:
// at a good prime P of F, Frob has char poly  x^2 - a_P x + Norm(P)  with
// a_P = a_P(f) in O_{K_f}.  The degree-4 (over Q) L-polynomial of B_f at P is
//   Norm_{K_f/Q}( 1 - a_P T + Norm(P) T^2 )   (RM => the two Q-conjugate factors).
// This is the exact Frobenius fingerprint any genus-2 model of B_f must reproduce.
//
// Run: magma goal2_mod5_Bf.m
SetColumns(0);
d := 2; F<w> := QuadraticField(d); OF := Integers(F);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);

LOGF := "goal2_mod5_Bf_out.txt";
PrintFile(LOGF, "# goal2_mod5_Bf: data of f and its abelian surface B_f" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1];
nlev := P2^3 * P5^2;
LOG(Sprintf("level 𝔫 = P2^3*P5^2, norm %o", Norm(nlev)));
t0 := Cputime();
M := HilbertCuspForms(F, nlev, [2,2]);
Dc := NewformDecomposition(NewSubspace(M));
LOG(Sprintf("space + decomposition built (%.1o s), %o orbits", Cputime(t0), #Dc));

ef := Eigenform(Dc[6]);
K := HeckeEigenvalueField(Dc[6]); OK := Integers(K);
LOG(Sprintf("orbit 6 Hecke field K_f: degree %o", Degree(K)));
LOG(Sprintf("  defining poly: %o", DefiningPolynomial(K)));
if Degree(K) eq 2 then
    D := Discriminant(Integers(K));
    LOG(Sprintf("  K_f = Q(sqrt %o), disc %o", Squarefree(Numerator(D)*Denominator(D)), D));
    LOG(Sprintf("  totally real: %o", IsTotallyReal(K)));
end if;

// Hecke eigenvalues at good primes of F (norm up to bound), and the B_f L-polys.
LOG("");
LOG("# P(norm) : a_P(f) in K_f : B_f L-poly over Q = Norm_{K/Q}(1 - a_P T + NP T^2)");
RT<T> := PolynomialRing(QQ);
badN := {2,5};
cnt := 0;
for l in PrimesInInterval(2, 60) do
    for pf in Factorization(l*OF) do
        P := pf[1]; NP := Norm(P);
        if l in badN then continue; end if;
        aP := HeckeEigenvalue(ef, P);
        // char poly of Frob on B_f over F is x^2 - aP x + NP (in K_f[x]);
        // its Q-L-polynomial is the K/Q-norm of (1 - aP T + NP T^2).
        if Degree(K) eq 2 then
            // Norm_{K/Q}(NP T^2 - aP T + 1) = NP^2 T^4 - NP*s T^3 + (2NP+nm) T^2 - s T + 1
            // with s = Tr_{K/Q}(aP), nm = Nm_{K/Q}(aP).
            s := Trace(K!aP); nm := Norm(K!aP);
            LqQ := NP^2*T^4 - NP*s*T^3 + (2*NP + nm)*T^2 - s*T + 1;
            LOG(Sprintf("P(%o) : %o : %o", NP, aP, LqQ));
        else
            LOG(Sprintf("P(%o) : %o", NP, aP));
        end if;
        cnt +:= 1;
    end for;
end for;
LOG(Sprintf("\n%o good primes tabulated. This is the Frobenius fingerprint of B_f.", cnt));
exit;
