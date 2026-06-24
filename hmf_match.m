// hmf_match.m -- Match the 2-dim mod-11 reps of Curve A against Hilbert modular forms via
// the (trace, norm) fingerprint, allowing quadratic Hecke fields.
//
// Target (from frobdata, every prime): with h_P = x^4 + c1 x^3 + c2 x^2 + c3 x + c4,
//     targetTr = -c1 mod 11   (= tr rho_1 + tr rho_2)
//     targetNm = c2 - 2 N(P) mod 11  (= tr rho_1 * tr rho_2)
// A Hilbert newform f over Q(sqrt5) (parallel wt 2, Hecke field F) matches if for every
// good K-prime P:  Tr_{F/Q}(a_P(f)) = targetTr  and  Nm_{F/Q}(a_P(f)) = targetNm  (mod 11).
// (When F is real-quadratic with 11 split, a_P(f) reduces to the pair {tr rho_1, tr rho_2}.)
//
// This re-does the small-level scan with the CORRECT target (the earlier hmf_smalllevel_scan
// used the conflated single trace). If a small (2-power) level matches, rho_i is unramified
// at 66179 (level-lowering) and the heavy level-66179 computation is avoided.
//
// Run: magma hmf_match.m

SetColumns(0);
load "mod11_frobdata.txt";
K<w> := QuadraticField(5); OK := Integers(K);
F11 := GF(11);
PBound := 200;

function PrimeFromTag(p, rid)
    if rid eq -1 then return ideal<OK | p>; end if;
    return ideal<OK | p, w - rid>;
end function;

// ---- build target fingerprint, write to file ----
target := [];   // <P-ideal, p, targetTr, targetNm>
trnm_rows := [];
for row in frobdata do
    Np := row[1]; p := row[2]; rid := row[3]; cs := row[4];
    tTr := F11 ! (-cs[2]);
    tNm := F11 ! (cs[3] - 2*Np);
    Append(~trnm_rows, <Np, p, rid, Integers()!tTr, Integers()!tNm>);
    if p le PBound and p ne 11 and p ne 2 then
        Append(~target, <PrimeFromTag(p,rid), p, tTr, tNm>);
    end if;
end for;
PrintFile("mod11_trnm.txt", "// Curve A: (trace,norm) fingerprint of {rho_1,rho_2}, mod 11" : Overwrite := true);
PrintFile("mod11_trnm.txt", "// rows <N(P), p, rid, Tr mod 11, Nm mod 11>");
PrintFile("mod11_trnm.txt", Sprintf("trnm := %o;", trnm_rows));
printf "Wrote mod11_trnm.txt (%o rows). Using %o K-primes (norm-p <= %o) for matching.\n\n",
    #trnm_rows, #target, PBound;

// ---- reduce a Hecke eigenvalue's Tr_{F/Q}, Nm_{F/Q} mod 11 ----
function TrNmMod11(a)
    F := Parent(a);
    if F cmpeq Rationals() or Type(F) eq FldRat then
        return F11!a, F11!a;            // degree-1: Tr=Nm=a (won't match generically)
    end if;
    return F11!(Integers()!Trace(a)), F11!(Integers()!Norm(a));
end function;

// ---- scan all {2,3,5}-supported levels (the level if rho_i lowers at 66179) ----
// A is bad only within {2,3,5,11,66179}; lowering at 66179 leaves level supported on
// {2,3,5}. 2 inert (N=4), 3 inert (N=9), 5 ramified (N=5).
p2 := Factorization(2*OK)[1][1];
p3 := Factorization(3*OK)[1][1];
p5 := Factorization(5*OK)[1][1];
NormCap := 6000;
levels := [];
for a in [0..7], b in [0..4], c in [0..4] do
    N := p2^a * p3^b * p5^c;
    if Norm(N) ge 31 and Norm(N) le NormCap then Append(~levels, N); end if;
end for;
Sort(~levels, func<x,y | Norm(x)-Norm(y)>);
printf "Scanning %o levels supported on {2,3,5}, norm in [31,%o]...\n\n", #levels, NormCap;

anyfull := false;
for N in levels do
    M := HilbertCuspForms(K, N, [2,2]);
    Mnew := NewSubspace(M);
    if Dimension(Mnew) eq 0 then continue; end if;
    decomp := NewformDecomposition(Mnew);
    bestorbit := 0;
    for nf in decomp do
        eig := Eigenform(nf);
        Fdeg := Degree(HeckeEigenvalueField(nf));
        matched := 0; okall := true;
        for tP in target do
            if Gcd(Norm(N), tP[2]) ne 1 then continue; end if;   // skip primes dividing level
            tr, nm := TrNmMod11(HeckeEigenvalue(eig, tP[1]));
            if tr ne tP[3] or nm ne tP[4] then okall := false; break; end if;
            matched +:= 1;
        end for;
        if matched gt bestorbit then bestorbit := matched; end if;
        if okall and matched ge 20 then
            printf "level norm %o, orbit dim %o (Hecke deg %o): *** FULL MATCH (%o primes) ***\n",
                Norm(N), Dimension(nf), Fdeg, matched;
            anyfull := true;
        end if;
    end for;
    printf "level norm %-5o: dim %o, %o orbits, best prefix-match %o/%o\n",
        Norm(N), Dimension(Mnew), #decomp, bestorbit, #target;
end for;
printf "\n%o\n", anyfull select "FOUND a small-level match -- level-lowering at 66179!"
    else "No {2,3,5}-level match. rho_i is genuinely ramified at 66179 (level >= 66179).";
