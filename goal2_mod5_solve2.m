// goal2_mod5_solve2.m -- robust knapsack CRT reconstruction of B_f's H_12 point (e,f).
// vs solve.m: (1) verify every candidate against ALL cached primes (score/#verify),
// allowing SLACK bad-reduction primes; (2) LEAVE-ONE-OUT over the CRT pool, to catch the
// case where the true point is bad at one CRT prime (its residue not in that M_l).
// Leave-one-out is cheap: dropping prime l costs base/|M_l| combos, so all drops together
// ~ base * sum(1/|M_l|) ~ 1.5x base.
//
// Args: BUDGET (combo cap for the FULL pool), SLACK (allowed bad primes, default 2),
//       LOO ("1" to also run leave-one-out, default "1").
// Run: magma BUDGET:=20000000 SLACK:=2 LOO:=1 goal2_mod5_solve2.m
SetColumns(0);
if assigned BUDGET then BUDGET := StringToInteger(BUDGET); else BUDGET := 20000000; end if;
if assigned SLACK then SLACK := StringToInteger(SLACK); else SLACK := 2; end if;
if assigned LOO then LOO := LOO eq "1"; else LOO := true; end if;

load "mcache.m";
LOGF := "goal2_mod5_solve2_out.txt";
PrintFile(LOGF, Sprintf("# robust knapsack CRT + LOO, BUDGET=%o SLACK=%o LOO=%o", BUDGET, SLACK, LOO) : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

all := [7,17,23,31,41,47,71,73,79,89,97,103,113,127,137,151,167,191,193,199];
all := [ l : l in all | IsDefined(Mcache, l) and #Mcache[l] gt 0 ];
sz := AssociativeArray(); for l in all do sz[l] := #Mcache[l]; end for;
MS := AssociativeArray();
for l in all do MS[l] := { <x[1],x[2]> : x in Mcache[l] }; end for;
function memOK(e0,f0,l)
    if (Denominator(e0) mod l eq 0) or (Denominator(f0) mod l eq 0) then return false; end if;
    return <Integers()!(GF(l)!e0), Integers()!(GF(l)!f0)> in MS[l];
end function;

// knapsack: greedy by log(l)/log(|M_l|) under BUDGET
byratio := Sort(all, func< a,b | (Log(b)/Log(sz[b])) - (Log(a)/Log(sz[a])) gt 0 select 1 else -1 >);
pool := []; prod := 1;
for l in byratio do
    if prod * sz[l] le BUDGET then Append(~pool, l); prod *:= sz[l]; end if;
end for;
pool := Sort(pool);
LOG(Sprintf("pool = %o  |M|=%o  combos=%o", pool, [sz[l]:l in pool], prod));

// scan one CRT subset S; verify against `all \ S`, keep score >= #verify - SLACK
allbest := 0; survivors := [];
procedure scanSubset(S, ~allbest, ~survivors)
    M := &*S; H := Isqrt(M div 2); ZM := Integers(M);
    verify := [ l : l in all | l notin S ];
    np := #S; Ms := [ Mcache[l] : l in S ]; sizes := [ #m : m in Ms ]; total := &*sizes;
    co := [ Integers() | ];
    for j in [1..np] do l := S[j]; Mj := M div l; Append(~co, Mj*(InverseMod(Mj mod l, l))); end for;
    LOG(Sprintf("  scan S=%o (combos %o, reach ~%o, verify %o primes)", S, total, H, #verify));
    t0 := Cputime(); tested := 0; bsc := 0;
    for cc in [0..total-1] do
        t := cc; Es := 0; Fs := 0;
        for j in [1..np] do
            idx := (t mod sizes[j]) + 1; t := t div sizes[j]; pr := Ms[j][idx];
            Es +:= pr[1]*co[j]; Fs +:= pr[2]*co[j];
        end for;
        be, qe := RationalReconstruction(ZM ! Es); if not be then continue; end if;
        bf, qf := RationalReconstruction(ZM ! Fs); if not bf then continue; end if;
        tested +:= 1;
        sc := 0; for l in verify do if memOK(qe,qf,l) then sc +:= 1; end if; end for;
        if sc gt bsc then bsc := sc; end if;
        if sc ge #verify - SLACK then
            Append(~survivors, <qe,qf,sc,#verify>);
            LOG(Sprintf("    *** CANDIDATE (e,f)=(%o,%o) score %o/%o [drop %o] ***", qe, qf, sc, #verify, [l:l in all|l notin S]));
        end if;
    end for;
    if bsc gt allbest then allbest := bsc; end if;
    LOG(Sprintf("    done: %o combos %.1o s, %o RR-valid, best score %o/%o", total, Cputime(t0), tested, bsc, #verify));
end procedure;

// full pool
scanSubset(pool, ~allbest, ~survivors);
// leave-one-out
if LOO then
    for i in [1..#pool] do
        S := [ pool[j] : j in [1..#pool] | j ne i ];
        scanSubset(S, ~allbest, ~survivors);
    end for;
end if;

LOG(Sprintf("\n=== GLOBAL best score %o ; %o survivor(s) ===", allbest, #survivors));
for s in survivors do LOG(Sprintf("SURVIVOR (e,f)=(%o,%o) score %o/%o", s[1],s[2],s[3],s[4])); end for;
exit;
