// identify_f.m -- identify the char-0 Hilbert newform f at level Pbad (norm 66179),
// parallel weight 2, whose mod-11 reduction is the eigenform with a_P^2 in {t1^2,t2^2}.
//
// We know f's mod-11 Hecke eigenvalues (refine_serre5_signs.txt). Decompose S_2(Pbad) into
// newform orbits, and find the orbit that reduces mod 11 (via some prime above 11 in its Hecke
// field) to our measured signature. Then report: Hecke field, CM, and small a_P.
//
// Run: magma identify_f.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);

LOGF := "identify_f_out.txt";
PrintFile(LOGF, "# identify_f" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

QI_A := [K| 1, 1/1459240*(243125*w - 482787),
    1/2229718720*(-54798934*w + 127710391),
    1/8517525510400*(182422196780*w - 406668692089),
    1/65073894899456000*(38134182372761*w - 85270624049895) ];
Pbad := [P : P in BadPrimesFromInvariants(QI_A, K) | Norm(P) eq 66179][1];
function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;

// measured mod-11 signature of f at a few primes <N,p,rid,a mod 11>
sig := [ <49,7,-1,7>, <29,29,18,7>, <29,29,11,5>, <31,31,6,7>, <31,31,25,1>,
         <41,41,28,0>, <41,41,13,4>, <59,59,8,3>, <59,59,51,7>, <71,71,54,3> ];
sigP := [ <PrimeFromTag(s[2],s[3]), F11!s[4]> : s in sig ];

LOG(Sprintf("level Pbad norm %o", Norm(Pbad)));
M := HilbertCuspForms(K, Pbad, [2,2]);
LOG(Sprintf("cuspidal dim %o; decomposing into newforms...", Dimension(M)));
t0 := Cputime();
NS := NewSubspace(M);
D := NewformDecomposition(NS);
LOG(Sprintf("  %o newform orbits, decomposition time %o", #D, Cputime(t0)));

// for each orbit, get eigenform, Hecke field, and reduce a_P mod (a prime above 11) to match sig
for idx in [1..#D] do
    f := Eigenform(D[idx]);
    // Hecke field
    a7 := HeckeEigenvalue(f, sigP[1][1]);
    E := Parent(a7);
    dE := (Type(E) eq FldRat) select 1 else Degree(E);
    // primes above 11 in E
    if dE eq 1 then
        lambdas := [* 11 *];   // reduce rationals mod 11
    else
        OE := Integers(E);
        lambdas := [* Pl[1] : Pl in Factorization(11*OE) *];
    end if;
    matched := false;
    for lam in lambdas do
        ok := true;
        for sp in sigP do
            ap := HeckeEigenvalue(f, sp[1]);
            if dE eq 1 then
                red := F11!(Integers()!ap);   // ap rational integer (level prime to 11)
            else
                Fq, redmap := ResidueClassField(lam);
                if #Fq ne 11 then ok := false; break; end if;   // need residue field F_11
                red := F11!(redmap(OE!ap));
            end if;
            if red ne sp[2] then ok := false; break; end if;
        end for;
        if ok then matched := true; break; end if;
    end for;
    if matched then
        LOG(Sprintf("*** orbit %o MATCHES f. Hecke field degree %o ***", idx, dE));
        if dE eq 1 then LOG("  Hecke field = Q (rational eigenform!)");
        else LOG(Sprintf("  Hecke field E = %o", DefiningPolynomial(E))); end if;
        // CM test: count a_P = 0 over many primes
        nz := 0; ntot := 0; bcdiff := 0; nsplit := 0;
        for N := 1 to 60 do end for;
        // sample primes
        prs := PrimesUpTo(120, K);
        for P in prs do
            if Norm(P + Pbad) ne 1 or Norm(P) eq 11 then continue; end if;
            if Norm(P) gt 200 then continue; end if;
            ap := HeckeEigenvalue(f, P);
            if ap eq 0 then nz +:= 1; end if;
            ntot +:= 1;
        end for;
        LOG(Sprintf("  a_P = 0 at %o / %o primes (N(P)<=200) -- CM would be ~50%%", nz, ntot));
        // base change test: compare a_P, a_Pbar for split rational primes
        for p in [29,31,41,59,71,89] do
            fac := [Pp[1] : Pp in Factorization(p*OK) | Norm(Pp[1]) eq p];
            if #fac eq 2 then
                a := HeckeEigenvalue(f, fac[1]); b := HeckeEigenvalue(f, fac[2]);
                nsplit +:= 1; if a ne b then bcdiff +:= 1; end if;
            end if;
        end for;
        LOG(Sprintf("  base-change-from-Q: a_P != a_Pbar at %o/%o split primes (>0 => NOT base change)", bcdiff, nsplit));
        // print a few char-0 eigenvalues
        for sp in sigP[1..6] do
            LOG(Sprintf("    a_P(N=%o) = %o", Norm(sp[1]), HeckeEigenvalue(f, sp[1])));
        end for;
    end if;
end for;
LOG("DONE.");
exit;
