// goal2_mod5_solve.m -- reconstruct B_f's H_12 moduli point (e,f) in Q^2 by CRT over a
// KNAPSACK-selected subset of split primes: maximize the modulus prod(l) (=> rational-
// reconstruction height reach ~ sqrt(M/2)) subject to a brute-force combo budget
// prod(|M_l|) <= BUDGET.  Candidates are then membership-verified at the primes NOT used
// in the CRT (the true point passes all with good reduction).
//
// Uses cached match-sets from mcache.m (built by goal2_mod5_mbuild.m).
// Optional arg: BUDGET (combo cap).  Run: magma BUDGET:=20000000 goal2_mod5_solve.m
SetColumns(0);
if assigned BUDGET then BUDGET := StringToInteger(BUDGET); else BUDGET := 20000000; end if;

load "mcache.m";   // Mcache[l] = [ <e,f>, ... ]  (integer reps in [0,l))

// fingerprint <s,n> for membership verification
target := AssociativeArray();
for t in [ <7,-2,-2>,<17,6,-3>,<23,14,46>,<31,6,6>,<41,-8,-11>,<47,0,-48>,<71,4,-8>,
           <73,-8,13>,<79,2,-2>,<89,-6,-3>,<97,-4,-8>,<103,16,16>,<113,32,253>,
           <127,18,-66>,<137,-26,121>,<151,14,-98>,<167,-12,-72>,<191,6,-138>,
           <193,34,181>,<199,8,-176> ] do target[t[1]]:=<t[2],t[3]>; end for;

LOGF := "goal2_mod5_solve_out.txt";
PrintFile(LOGF, Sprintf("# knapsack CRT solve, BUDGET=%o", BUDGET) : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

all := [7,17,23,31,41,47,71,73,79,89,97,103,113,127,137,151,167,191,193,199];
all := [ l : l in all | IsDefined(Mcache, l) and #Mcache[l] gt 0 ];
sz := AssociativeArray(); for l in all do sz[l] := #Mcache[l]; end for;

// membership set (as a Magma set) per prime
MS := AssociativeArray();
for l in all do MS[l] := { <x[1],x[2]> : x in Mcache[l] }; end for;
function memOK(e0,f0,l)
    if (Denominator(e0) mod l eq 0) or (Denominator(f0) mod l eq 0) then return false; end if;
    return <Integers()!(GF(l)!e0), Integers()!(GF(l)!f0)> in MS[l];
end function;

// --- knapsack: greedy by log(l)/log(|M_l|), fill while prod(|M|) <= BUDGET ---
byratio := Sort(all, func< a,b | (Log(b)/Log(sz[b])) - (Log(a)/Log(sz[a])) gt 0 select 1 else -1 >);
S := []; prod := 1;
for l in byratio do
    if sz[l] eq 1 then Append(~S, l); continue; end if;
    if prod * sz[l] le BUDGET then Append(~S, l); prod *:= sz[l]; end if;
end for;
S := Sort(S);
fresh := [ l : l in all | l notin S ];
M := &*S; H := Isqrt(M div 2);
LOG(Sprintf("CRT primes S = %o", S));
LOG(Sprintf("  |M_l| = %o", [sz[l] : l in S]));
LOG(Sprintf("  combos prod|M| = %o", prod));
LOG(Sprintf("  modulus M = %o  (~%o digits) => height reach ~ %o", M, #IntegerToString(M), H));
LOG(Sprintf("verify primes (not in CRT) = %o", fresh));

ZM := Integers(M);
Ms := [ Mcache[l] : l in S ]; np := #S; sizes := [ #m : m in Ms ];
total := &*sizes;

// precompute per-prime CRT basis to speed inner loop:  x = sum_j r_j * co_j mod M
// where co_j = (M/l_j) * ((M/l_j)^{-1} mod l_j)
co := [ Integers() | ];
for j in [1..np] do
    l := S[j]; Mj := M div l;
    Append(~co, Mj * (InverseMod(Mj mod l, l)));
end for;

LOG("scanning...");
found := []; t0 := Cputime(); tested := 0; bestsc := 0;
for cc in [0..total-1] do
    t := cc; Esum := 0; Fsum := 0;
    for j in [1..np] do
        idx := (t mod sizes[j]) + 1; t := t div sizes[j];
        pr := Ms[j][idx];
        Esum +:= pr[1]*co[j];
        Fsum +:= pr[2]*co[j];
    end for;
    be, qe := RationalReconstruction(ZM ! Esum);
    if not be then continue; end if;
    bf, qf := RationalReconstruction(ZM ! Fsum);
    if not bf then continue; end if;
    tested +:= 1;
    // full fresh-score; allow a few bad-reduction fresh primes.
    sc := 0;
    for l in fresh do if memOK(qe, qf, l) then sc +:= 1; end if; end for;
    if sc gt bestsc then bestsc := sc; end if;
    if sc ge #fresh - 2 then
        Append(~found, <qe,qf,sc>);
        LOG(Sprintf("  *** CANDIDATE (e,f) = (%o, %o)  fresh-score %o/%o ***", qe, qf, sc, #fresh));
    end if;
    if cc mod 2000000 eq 0 and cc gt 0 then
        LOG(Sprintf("  ...%o / %o combos [%.1o s], %o RR-valid, best fresh-score %o/%o",
                    cc, total, Cputime(t0), tested, bestsc, #fresh));
    end if;
end for;
LOG(Sprintf("scan done: %o combos in %.1o s, %o RR-valid candidates, best fresh-score %o/%o, %o survivors(>=%o)",
            total, Cputime(t0), tested, bestsc, #fresh, #found, #fresh-2));
for ef in found do LOG(Sprintf("SURVIVOR (e,f) = (%o, %o)  score %o", ef[1], ef[2], ef[3])); end for;
exit;
