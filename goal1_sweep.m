// goal1_sweep.m -- Goal 1: for each entry in sorted_output.jsonl, find a
// candidate Hilbert newform f over Q(sqrt d) whose induction matches A[p].
//
// Method (validated on the mod-7 Q(sqrt10) example, data idx 33):
//   A[p] is induced from F; at a rational prime l SPLIT in F (l = Pa*Pb),
//   trace of Frob_l on Ind(rho_f) = a_Pa(f) + a_Pb(f), which must equal
//   a_l(A) mod p.  (At inert primes the induced trace is 0 automatically.)
//   For each candidate level ideal n (Norm n = serre_conductor), decompose
//   the new cusp space into orbits and, reducing each orbit mod every prime
//   lambda|p of its Hecke field, test the split-prime trace identity.
//
// Resume-aware: entries already recorded in RESF are skipped.
// Tune CUTOFF to control how large a level we attempt.
//
// Run: magma goal1_sweep.m

SetColumns(0);
load "goal1_data.m";           // defines `entries`
LVLCUTOFF := 12000;           // skip entries whose FINAL level norm exceeds this
MAXIDEAL := 8;                 // cap candidate level ideals tried per entry
NSPLIT   := 20;                // split test primes to use
PRIMEBD  := 400;               // search split primes up to here
RESF     := "goal1_results.txt";

// ---- resume: collect already-done indices ----
done := {};
if OpenTest(RESF, "r") then
    fh := Open(RESF, "r");
    while true do
        s := Gets(fh);
        if IsEof(s) then break; end if;
        if #s ge 5 and s[1..5] eq "ENTRY" and Position(s, "DEFERRED") eq 0 then
            t := Split(s, " \t");
            done := done join { StringToInteger(t[2]) };
        end if;
    end while;
    delete fh;
end if;
if not OpenTest(RESF, "r") then
    PrintFile(RESF, "# Goal-1 sweep results: ENTRY <idx> d<D> p<P> serre<N> -> <status>");
end if;
procedure REC(s) printf "%o\n", s; PrintFile(RESF, s); end procedure;

// ---- ideals of a given (rational) norm ----
function ExpSolutions(fs, e)
    if #fs eq 0 then return (e eq 0) select [[Integers()|]] else [Parent([Integers()|])|]; end if;
    f1 := fs[1]; rest := fs[2..#fs]; sols := []; a := 0;
    while a*f1 le e do
        for tail in ExpSolutions(rest, e - a*f1) do Append(~sols, [a] cat tail); end for;
        a +:= 1;
    end while;
    return sols;
end function;

function IdealsOfNorm(OF, N)
    result := [ ideal<OF|1> ];
    for pe in Factorization(N) do
        ell := pe[1]; e := pe[2];
        Ps := [ pl[1] : pl in Factorization(ell*OF) ];
        fs := [ Valuation(Norm(P), ell) : P in Ps ];
        locals := [];
        for sol in ExpSolutions(fs, e) do
            J := ideal<OF|1>;
            for i in [1..#Ps] do J := J * Ps[i]^sol[i]; end for;
            Append(~locals, J);
        end for;
        result := [ I*J : I in result, J in locals ];
    end for;
    return result;
end function;

// ---- main loop ----
for ent in entries do
    idx := ent[1]; d := ent[2]; p := ent[3]; serre := ent[4];
    fc := ent[5]; hc := ent[6]; cond := ent[7];
    hdr := Sprintf("ENTRY %o d%o p%o serre%o", idx, d, p, serre);
    if idx in done then continue; end if;

    F<w> := QuadraticField(d); OF := Integers(F);
    QQ := Rationals(); Px<x> := PolynomialRing(QQ);
    fpol := &+[ fc[i]*x^(i-1) : i in [1..#fc] ];
    hpol := &+[ hc[i]*x^(i-1) : i in [1..#hc] ];
    C := HyperellipticCurve(fpol, hpol);
    dC := Integers()!Discriminant(C);
    Fp := GF(p);

    // Level construction. Serre conductor is prime-to-p.  If p is a BAD prime
    // (p | cond) it must be restored into the level: empirically the conductor
    // exponent at the prime above p is 2 (validated on idx 0,5,8 for p=5 inert).
    // If p is good (e.g. the mod-7 Q(sqrt10) family), level = SerreIdeal.
    pbad := (cond mod p eq 0);
    cids := IdealsOfNorm(OF, serre);
    if pbad then
        Pp := Factorization(p*OF)[1][1];
        cids := [ nn * Pp^2 : nn in cids ];
    end if;
    cids := [ nn : nn in cids | Norm(nn) le LVLCUTOFF ];
    if #cids eq 0 then
        REC(Sprintf("%o -> DEFERRED (final level norm > %o; pbad=%o)", hdr, LVLCUTOFF, pbad));
        continue;
    end if;

    // A fingerprint at split good primes: <l, Pa, Pb, a_l(A) mod p>
    targ := [];
    for l in PrimesInInterval(3, PRIMEBD) do
        if #targ ge NSPLIT then break; end if;
        if l eq p or (d mod l eq 0) or (dC mod l eq 0) then continue; end if;
        if not IsSquare(GF(l)!d) then continue; end if;       // split test
        fp := Factorization(l*OF);
        Cl := ChangeRing(C, GF(l));
        al := l + 1 - #Cl;
        Append(~targ, < l, fp[1][1], fp[2][1], Integers()!(Fp!al) >);
    end for;

    if #cids gt MAXIDEAL then cids := cids[1..MAXIDEAL]; end if;

    matched := false;
    t0 := Cputime();
    for ni -> nn in cids do
        M := HilbertCuspForms(F, nn, [2,2]);
        D := NewformDecomposition(NewSubspace(M));
        for oi in [1..#D] do
            e := Eigenform(D[oi]);
            // gather (a_Pa, a_Pb) in Hecke field for usable split primes
            E := HeckeEigenvalueField(D[oi]);
            isQ := Type(E) eq FldRat;
            if isQ then lams := [* <1> *]; else
                OE := Integers(E); lams := [* <pl[1]> : pl in Factorization(p*OE) *];
            end if;
            pairs := [];
            for t in targ do
                if Norm(nn) mod t[1] eq 0 then continue; end if;
                Append(~pairs, < HeckeEigenvalue(e,t[2]), HeckeEigenvalue(e,t[3]), t[4] >);
            end for;
            if #pairs lt 6 then continue; end if;
            for lam in lams do
                if isQ then Fq := Fp; red := func< a | Fp!(Integers()!a) >;
                else Fq, r0 := ResidueClassField(lam[1]); red := func< a | r0(Integers(E)!a) >; end if;
                nbad := 0;
                for pr in pairs do
                    if red(pr[1]) + red(pr[2]) ne Fq!(pr[3]) then nbad +:= 1; end if;
                end for;
                if nbad eq 0 then
                    dEg := isQ select 1 else Degree(E);
                    REC(Sprintf("%o -> MATCH  ideal#%o(norm %o)  orbit %o  HeckeDeg %o  res F_%o  tests %o/%o  [%.1os]",
                        hdr, ni, Norm(nn), oi, dEg, #Fq, #pairs, #pairs, Cputime(t0)));
                    matched := true; break;
                end if;
            end for;
            if matched then break; end if;
        end for;
        if matched then break; end if;
    end for;
    if not matched then
        REC(Sprintf("%o -> NO MATCH  (%o ideals, %.1os)", hdr, #cids, Cputime(t0)));
    end if;
end for;
REC("# pass done (CUTOFF " cat IntegerToString(CUTOFF) cat ")");
exit;
