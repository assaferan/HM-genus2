// identify_g.m -- identify the level-Pbad newform g with rho_g = rho_i (x) chi, chi = chi_g*chi_{-1}.
// Decompose S_2(Pbad) into newform orbits; for each, reduce its Hecke eigenvalues mod a prime
// above 11 and test whether a_P is a root of the chi-twisted fingerprint y^2 - chi(P)Tr y + Nm.
//
// Run: magma identify_g.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;
import "Reduction.m": MinimalTwist, IntegralModelOf, BadPrimesFromInvariants;
K<s5> := QuadraticField(5); OK := Integers(K); PK<x> := PolynomialRing(K); F11 := GF(11);
LOGF := "identify_g_out.txt"; PrintFile(LOGF, "# identify_g" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI := [K| 1, 1/1459240*(243125*s5 - 482787), 1/2229718720*(-54798934*s5 + 127710391),
    1/8517525510400*(182422196780*s5 - 406668692089), 1/65073894899456000*(38134182372761*s5 - 85270624049895) ];
C := Genus2CurveFromIgusa(QI, K); Cmin := MinimalTwist(C, K); Cint := IntegralModelOf(Cmin);
ff,hh := HyperellipticPolynomials(Cint);
F0 := PK![OK!c : c in Coefficients(ff + (Parent(ff)!hh)^2/4)];
Msq := &*([ideal<OK|1>] cat [pe[1]^(pe[2] div 2) : pe in Factorization(&+[ideal<OK|c> : c in Coefficients(F0)|c ne 0])]);
_, alpha := IsPrincipal(Msq);
F1 := PK![ (OK!c)/alpha^2 : c in Coefficients(F0) ];
_, g := IsPrincipal(&+[ideal<OK|c> : c in Coefficients(F1)|c ne 0]);
Pbad := [P : P in BadPrimesFromInvariants(QI, K) | Norm(P) eq 66179][1];
function PrimeFromTag(p, rid) if rid eq -1 then return ideal<OK|p>; end if; return ideal<OK|p, s5-rid>; end function;
function ChiG(P) Fp,red := ResidueClassField(P); v := red(OK!g); return (v ne 0 and IsSquare(v)) select 1 else -1; end function;
function ChiM1(P) Fp,red := ResidueClassField(P); return IsSquare(red(OK!(-1))) select 1 else -1; end function;  // chi_{-1}

load "mod11_trnm.txt";
targs := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    P := PrimeFromTag(row[2], row[3]);
    if Norm(P+Pbad) ne 1 then continue; end if;
    chi := ChiG(P)*ChiM1(P);   // combined twist sign in {+-1}
    Append(~targs, <P, chi, Integers()!(F11!row[4]), Integers()!(F11!row[5]), row[1]>);
end for;
Sort(~targs, func<a,b|a[5]-b[5]>);
targs := targs[1..14];   // 14 test primes

LOG(Sprintf("level Pbad norm %o; decomposing S_2 into newforms...", Norm(Pbad)));
M := HilbertCuspForms(K, Pbad, [2,2]);
t0 := Cputime(); D := NewformDecomposition(NewSubspace(M));
LOG(Sprintf("  %o orbits (%o s)", #D, Cputime(t0)));

for idx in [1..#D] do
    e := Eigenform(D[idx]);
    a := HeckeEigenvalue(e, targs[1][1]); E := Parent(a);
    dE := (Type(E) eq FldRat) select 1 else Degree(E);
    LOG(Sprintf("orbit %o: Hecke field degree %o", idx, dE));
    // primes above 11
    if dE eq 1 then
        lams := [* 0 *];   // sentinel for the rational prime 11
    else
        OE := Integers(E); lams := [* pl[1] : pl in Factorization(11*OE) *];
    end if;
    for li in [1..#lams] do
        lam := lams[li];
        if dE eq 1 then Fq := F11; fdeg := 1; else Fq, red := ResidueClassField(lam); fdeg := Degree(Fq); end if;
        ok := true; nbad := 0;
        for tg in targs do
            aP := HeckeEigenvalue(e, tg[1]);
            if dE eq 1 then aPb := Fq!(Integers()!aP); else aPb := red(Integers(E)!aP); end if;
            twTr := Fq!(tg[2] mod 11) * Fq!(tg[3]);   // chi(P)*Tr  (chi=+-1)
            val := aPb^2 - twTr*aPb + Fq!(tg[4]);
            if val ne 0 then ok := false; nbad +:= 1; end if;
        end for;
        if ok then
            LOG(Sprintf("  *** MATCH via prime above 11 with residue F_(11^%o): orbit %o IS g ***", fdeg, idx));
            if dE gt 1 then LOG(Sprintf("     Hecke field: %o", DefiningPolynomial(E)));
                LOG(Sprintf("     11 factors in E as: %o", [<Degree(ResidueClassField(pl[1])), pl[2]> : pl in Factorization(11*Integers(E))])); end if;
            LOG(Sprintf("     small eigenvalues a_P: %o", [<Norm(tg[1]), HeckeEigenvalue(e, tg[1])> : tg in targs[1..5]]));
        else
            LOG(Sprintf("  prime %o (F_11^%o): no match (%o/%o primes fail)", li, fdeg, nbad, #targs));
        end if;
    end for;
end for;
LOG("DONE."); exit;
