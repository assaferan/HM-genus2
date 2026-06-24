// build_match_hmf.m -- Find the Hilbert modular form(s) over Q(sqrt5) realizing the 2-dim
// mod-11 reps rho_1, rho_2 of Curve A. THIS IS THE SERRE-MODULARITY MATCH.
//
// For each Hilbert newform f of parallel weight 2 at the level (Curve A's bad prime above
// 66179, possibly times a 2-power) and each degree-1 prime lambda | 11 in f's Hecke field,
// set t_P = a_P(f) mod lambda in F_11.  Then rho-bar_{f,lambda} = rho_1 (or rho_2) iff
//     t_P^2 - targetTr_P * t_P + targetNm_P = 0  (mod 11)   for every test prime P,
// where (targetTr, targetNm) = (-A, B-2N(P)) is the fingerprint in mod11_trnm.txt.
// Finding such an f IS evidence for Serre's conjecture over Q(sqrt5) (which is open there).
//
// Run: magma build_match_hmf.m

SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
import "Reduction.m": BadPrimesFromInvariants;

K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
load "mod11_trnm.txt";   // trnm := [ <Np,p,rid,Tr,Nm>, ... ]

plog := "match_hmf_progress.txt";
PrintFile(plog, "# HMF match progress" : Overwrite := true);
procedure LOG(s) PrintFile(plog, s); end procedure;

// Curve A invariants -> its bad prime above 66179
QI_A := [K|
    1,
    1/1459240*(243125*w - 482787),
    1/2229718720*(-54798934*w + 127710391),
    1/8517525510400*(182422196780*w - 406668692089),
    1/65073894899456000*(38134182372761*w - 85270624049895)
];
bad := BadPrimesFromInvariants(QI_A, K);
Pbad := [P : P in bad | Norm(P) eq 66179][1];
p2 := Factorization(2*OK)[1][1];
printf "Curve A bad prime: norm %o\n", Norm(Pbad);

function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;

// build test targets (skip 11)
targets := [];
for row in trnm do
    if row[2] eq 11 then continue; end if;
    Append(~targets, <PrimeFromTag(row[2],row[3]), F11!row[4], F11!row[5], row[2]>);
end for;
printf "%o test primes loaded.\n\n", #targets;

procedure CheckLevel(Level)
    LOG(Sprintf("level norm %o: building space...", Norm(Level)));
    M := HilbertCuspForms(K, Level, [2,2]);
    Mnew := NewSubspace(M);
    printf "=== level norm %o: dim_new = %o ===\n", Norm(Level), Dimension(Mnew);
    LOG(Sprintf("  dim %o; decomposing...", Dimension(Mnew)));
    D := NewformDecomposition(Mnew);
    printf "  %o newform orbits\n", #D;
    LOG(Sprintf("  %o orbits; checking each...", #D));
    found := 0;
    for idx in [1..#D] do
        nf := D[idx];
        eig := Eigenform(nf);
        F := HeckeEigenvalueField(nf);
        OF := Type(F) eq FldRat select Integers() else Integers(F);
        // degree-1 primes above 11 in F
        if Type(F) eq FldRat then
            lams := [11];   // F = Q: residue field F_11
        else
            lams := [fa[1] : fa in Factorization(11*OF) | #ResidueClassField(fa[1]) eq 11];
        end if;
        for lam in lams do
            if Type(F) eq FldRat then
                redmap := func<a | F11 ! (Integers()!a)>;
            else
                Ff, rr := ResidueClassField(lam);
                redmap := func<a | F11 ! (Integers()!rr(OF!a))>;
            end if;
            okall := true; nchk := 0;
            for tg in targets do
                if Norm(tg[1] + Level) ne 1 then continue; end if;   // coprime to level
                t := redmap(HeckeEigenvalue(eig, tg[1]));
                if t^2 - tg[2]*t + tg[3] ne 0 then okall := false; break; end if;
                nchk +:= 1;
            end for;
            if okall and nchk ge 20 then
                found +:= 1;
                printf "  *** MATCH: orbit %o (dim %o, Hecke deg %o), lambda|11 deg 1; verified %o primes ***\n",
                    idx, Dimension(nf), Degree(F), nchk;
                LOG(Sprintf("  *** MATCH orbit %o dim %o Hecke-deg %o (%o primes) ***", idx, Dimension(nf), Degree(F), nchk));
            end if;
        end for;
        if idx mod 20 eq 0 then LOG(Sprintf("  ...checked %o/%o orbits, %o matches", idx, #D, found)); end if;
    end for;
    printf "  => %o matching (newform, lambda) found at level norm %o\n\n", found, Norm(Level);
    LOG(Sprintf("  done level norm %o: %o matches", Norm(Level), found));
end procedure;

// level = bad prime (2-part trivial); if no match, try times 2-power.
for Lev in [Pbad, Pbad*p2, Pbad*p2^2] do
    CheckLevel(Lev);
end for;
