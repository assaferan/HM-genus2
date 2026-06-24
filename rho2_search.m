// rho2_search.m -- search for the second piece rho_2's Hilbert modular form at level p*(2)^2
// (norm 1058864, dim ~22059). Twist-robust V built via ker(T_P - lambda) nullspaces (no T^2).
// sigma_1 appears here as oldforms (mult 3); a NEW eigensystem with different a_P = rho_2.
//
// Designed for lovelace (256 cores, ~2TB). Logs progress to rho2_progress.txt.
//
// Run: magma rho2_search.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K); F11 := GF(11); R<y> := PolynomialRing(F11);
load "mod11_trnm.txt";

plog := "rho2_progress.txt";
PrintFile(plog, "# rho2 search at level p*(2)^2" : Overwrite := true);
procedure LOG(s) PrintFile(plog, s); printf "%o\n", s; end procedure;

QI_A := [K| 1, 1/1459240*(243125*w - 482787),
    1/2229718720*(-54798934*w + 127710391),
    1/8517525510400*(182422196780*w - 406668692089),
    1/65073894899456000*(38134182372761*w - 85270624049895) ];
Pbad := [P : P in BadPrimesFromInvariants(QI_A, K) | Norm(P) eq 66179][1];
p2 := Factorization(2*OK)[1][1];
Level := Pbad*p2^2;
function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if; return ideal<OK | p, w - rid>; end function;
targets := [];
for r in trnm do
    if r[2] eq 11 then continue; end if;
    Append(~targets, <PrimeFromTag(r[2],r[3]), F11!r[4], F11!r[5], r[1], r[2], r[3]>);
end for;
Sort(~targets, func<a,b | a[4]-b[4]>);

t0 := Cputime();
M := HilbertCuspForms(K, Level, [2,2]);
d := Dimension(M);
LOG(Sprintf("level norm %o, cuspidal dim = %o  (setup %os)", Norm(Level), d, Round(Cputime(t0))));

Vd := VectorSpace(F11, d);
B := 0;   // basis matrix (k x d) of the current candidate space, once initialized
built := 0;
for tg in targets do
    if Norm(tg[1] + Level) ne 1 then continue; end if;
    tp := Cputime();
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
    roots := { r[1] : r in Roots(y^2 - tg[2]*y + tg[3]) } join { r[1] : r in Roots(y^2 + tg[2]*y + tg[3]) };
    if Type(B) eq RngIntElt then
        // FIRST prime: full dim-d nullspaces (expensive, ~once)
        S := sub<Vd | >;
        for lam in roots do S := S + Kernel(T - lam*IdentityMatrix(F11, d)); end for;
        B := BasisMatrix(S);
    else
        // SUBSEQUENT primes (cheap): v=c*B is in ker(T-lam) iff c*(B*T - lam*B)=0.
        k := Nrows(B);
        BT := B*T;                               // k x d
        Vk := VectorSpace(F11, k);
        NC := sub<Vk | >;
        for lam in roots do NC := NC + Kernel(BT - lam*B); end for;   // left kernel, k-dim
        B := BasisMatrix(NC) * B;                // lift coords back to d-dim
    end if;
    T := 0;
    built +:= 1;
    LOG(Sprintf("  after N=%o (p=%o): dim V = %o  [%o roots, %os]",
        tg[4], tg[5], Nrows(B), #roots, Round(Cputime(tp))));
    if Nrows(B) eq 0 then LOG("  -> collapsed; rho_2 NOT here."); exit; end if;
    if built ge 14 then break; end if;
end for;
LOG(Sprintf("\nFINAL dim V = %o", Nrows(B)));

// extract eigensystems on V and compare a_P to sigma_1 (oldforms, mult 3) vs a NEW system (rho_2).
LOG("\nN(P)  p   {a_P on V}        {q-roots}");
cnt := 0;
for tg in targets do
    if Norm(tg[1]+Level) ne 1 then continue; end if;
    T := ChangeRing(Matrix(HeckeOperator(M, tg[1])), F11);
    fl, X := IsConsistent(B, B*T); T := 0;
    if not fl then LOG(Sprintf("%o %o  V not stable", tg[4], tg[5])); continue; end if;
    ev := { e[1] : e in Eigenvalues(X) };
    qr := { r[1] : r in Roots(y^2 - tg[2]*y + tg[3]) };
    LOG(Sprintf("%o %o  %o   %o", tg[4], tg[5], ev, qr));
    cnt +:= 1; if cnt ge 14 then break; end if;
end for;
LOG("DONE");
