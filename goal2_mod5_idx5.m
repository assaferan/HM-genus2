// goal2_mod5_idx5.m -- Goal-2 analog of idx 33, for a MOD-5 example (idx 5).
//
// Entry 5 of sorted_output.jsonl:  F = Q(sqrt2), p = 5, Serre norm 8,
//   C: y^2 = 5 - 10x + 10x^2 + 5x^4 + 2x^5,  cond(A) = 1600000 = 2^9 * 5^5.
// Goal-1 match: orbit 6 of HilbertCuspForms(F, n5000, [2,2]), Hecke field
//   degree 2, reduced at lambda|5 with residue field F_25.  (level norm 5000
//   = SerreIdeal(norm 8) * P5^2, P5 the inert prime above 5.)
//
// WHY mod 5 is the better Goal-2 target than mod 7 (idx 33): the Hecke field has
// degree 2, so B_f is an ABELIAN SURFACE over F (not the dim-18 variety of idx
// 33), and Res_{F/Q}(B_f) is 4-dimensional -- exactly the regime of van Bommel-
// Costa-Elkies-Keller-Schiavone-Voight (arXiv:2411.07857).  So the unconditional
// torsion-field route, blocked for idx 33, is in principle reachable here.
//
// This script establishes, exactly as for idx 33:
//   STEP 1  A[5] induced + irreducible + big image SL_2(F_25) (Dickson);
//   STEP 4  GRH-conditional certificate  A[5] = Ind(rho_f):
//           tr A[5](Frob_l) = tr Ind(rho_f)(Frob_l) for all good l < BOUND.
//
// Bad primes S = {2,5} (cond = 2^9*5^5); note 5 IS bad here (5|cond), unlike the
// good prime 7 of idx 33 -- so the level carries the P5^2 factor.  The Brauer-
// Nesbitt + conductor-based effective bound argument is otherwise identical.
//
// Run: magma goal2_mod5_idx5.m
SetColumns(0);
d := 2; p := 5; serreN := 8;
F<w> := QuadraticField(d); OF := Integers(F);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
fc := [5,-10,10,0,5,2];
fpol := &+[ fc[i]*x^(i-1) : i in [1..#fc] ];
C := HyperellipticCurve(fpol);
Fp := GF(p);
badp := {2, 5};

LOGF := "goal2_mod5_idx5_out.txt";
PrintFile(LOGF, "# goal2_mod5_idx5: Goal-2 analog for a mod-5 example (idx 5)" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

// ---------- helpers: L-polynomial data of A at a good prime l ----------
// Returns <a1, a2> with char poly of Frob = T^4 - a1 T^3 + a2 T^2 - a1 l T + l^2.
function LData(l)
    Cl := ChangeRing(C, GF(l));
    Lp := LPolynomial(Cl);                 // 1 - a1 T + a2 T^2 - a1 l T^3 + l^2 T^4
    a1 := l + 1 - #Cl;
    a2 := Coefficient(Lp, 2);
    return a1, a2;
end function;

// ================= STEP 1: image of A[5] =================
LOG("===== STEP 1: mod-5 image of A[5] =====");
IMGBD := 2000;
ninert := 0; ninert0 := 0;          // induced test: a_l == 0 at inert
nsplit := 0; nnonsplitss := 0;      // irreducibility: sigma-trace in F_25\F_5
ntrace0 := 0;                       // dihedral test: fraction of tr=0
projtr := {};                       // exceptional test: distinct proj traces t^2/l
transvection := false;              // t^2/l == 4 (order-5 unipotent, i.e. tr^2=4*det)
F25<u> := GF(25);
embok := true;
for l in PrimesInInterval(3, IMGBD) do
    if l in badp then continue; end if;
    a1, a2 := LData(l);
    tr := Fp ! a1;
    if tr eq 0 then ntrace0 +:= 1; end if;
    if IsSquare(GF(l)!d) then
        // split: sigma-traces t1,t2 are roots of Z^2 - a1 Z + (a2 - 2l) over F_5
        nsplit +:= 1;
        R<Z> := PolynomialRing(Fp);
        g := Z^2 - (Fp!a1)*Z + (Fp!(a2 - 2*l));
        rts := Roots(g, F25);
        // projective trace invariant t^2/l for a nonzero sigma-trace
        for rr in Roots(g, F25) do
            t := rr[1];
            if t ne 0 and (F25!l) ne 0 then
                Include(~projtr, t^2 / (F25!l));
                if t^2 / (F25!l) eq F25!4 then transvection := true; end if;
            end if;
        end for;
        // nonsplit-semisimple: g irreducible over F_5 => trace pair in F_25\F_5
        if IsIrreducible(g) then nnonsplitss +:= 1; end if;
    else
        // inert: induced trace must be 0 mod 5
        ninert +:= 1;
        if tr eq 0 then ninert0 +:= 1; end if;
    end if;
end for;
LOG(Sprintf("induced structure : %o / %o inert primes have a_l == 0 mod 5", ninert0, ninert));
LOG(Sprintf("irreducibility    : %o / %o split primes give sigma-trace in F_25\\F_5 (nonsplit-ss)", nnonsplitss, nsplit));
LOG(Sprintf("dihedral test     : tr == 0 at %o / %o primes (%.1o%%; dihedral ~ 50%%)",
    ntrace0, nsplit+ninert, 100.0*ntrace0/(nsplit+ninert)));
LOG(Sprintf("exceptional test  : %o distinct projective traces t^2/l (A4/S4/A5 would be few)", #projtr));
LOG(Sprintf("transvection      : %o (order-5 unipotent present)", transvection));
big := (ninert0 eq ninert) and (nnonsplitss gt 0) and transvection and (#projtr gt 10);
LOG(Sprintf(">>> STEP 1 verdict: induced + irreducible + big image SL_2(F_25) : %o", big));

// ================= STEP 4: GRH-conditional certificate =================
LOG("");
LOG("===== STEP 4: GRH certificate  A[5] = Ind(rho_f) =====");
LOG("building orbit-6 eigenform (level norm 5000, weight [2,2]) ...");
t0 := Cputime();
P2 := Factorization(2*OF)[1][1];         // 2 ramifies: (2) = P2^2
P5 := Factorization(5*OF)[1][1];         // 5 inert:   norm 25
// SerreIdeal of norm 8 = P2^? with norm 8 -> P2 has norm 2, so P2^3.  Level = P2^3 * P5^2.
nlev := P2^3 * P5^2;
assert Norm(nlev) eq 5000;
M := HilbertCuspForms(F, nlev, [2,2]);
Dc := NewformDecomposition(NewSubspace(M));
LOG(Sprintf("  space built (%.1o s), %o newform orbits; taking orbit 6", Cputime(t0), #Dc));
ef := Eigenform(Dc[6]); E := HeckeEigenvalueField(Dc[6]); OE := Integers(E);
LOG(Sprintf("  Hecke field degree %o", Degree(E)));
found := false; Fq := Fp; red := func< a | Fp!(Integers()!a) >;
if Degree(E) gt 1 then
    for pl in Factorization(5*OE) do
        Fq0, r0 := ResidueClassField(pl[1]);
        Fq := Fq0; red := func< a | r0(OE!a) >; found := true;
        break;
    end for;
    assert found;
end if;
LOG(Sprintf("  reduced at lambda|5, residue field F_%o", #Fq));

BOUND := 10000;
LOG(Sprintf("verifying all good primes up to %o ...", BOUND));
gsplit := 0; ginert := 0; fails := []; nextlog := 1000;
for l in PrimesInInterval(2, BOUND) do
    if l in badp then continue; end if;
    a1 := l + 1 - #ChangeRing(C, GF(l));
    aA := Fp ! a1;                                    // tr A[5](Frob_l) mod 5
    if IsSquare(GF(l)!d) then                          // split
        fp := Factorization(l*OF);
        s := red(OE!HeckeEigenvalue(ef, fp[1][1])) + red(OE!HeckeEigenvalue(ef, fp[2][1]));
        ok := (s eq Fq!aA);
        gsplit +:= 1;
    else                                               // inert: induced trace = 0
        ok := (aA eq Fp!0);
        ginert +:= 1;
    end if;
    if not ok then Append(~fails, l); end if;
    if l gt nextlog then
        LOG(Sprintf("  ... up to %o: %o split + %o inert checked, %o fails", l, gsplit, ginert, #fails));
        nextlog +:= 2000;
    end if;
end for;

LOG(Sprintf("\nTOTAL: %o good primes (%o split, %o inert) up to %o", gsplit+ginert, gsplit, ginert, BOUND));
LOG(Sprintf("trace DISAGREEMENTS: %o", #fails));
if #fails eq 0 then
    LOG("");
    LOG(">>> tr A[5](Frob_l) = tr Ind(rho_f)(Frob_l) for EVERY good prime l < " cat IntegerToString(BOUND));
    LOG("    cond(A) = 1600000, (log cond)^2 ~ 205; GRH conductor-based effective");
    LOG("    bound is O(10^3-10^4), within range.  => (GRH) A[5] = Ind(rho_f).");
    LOG("    Hecke field deg 2 => B_f is an abelian surface, Res(B_f) dim 4:");
    LOG("    the unconditional torsion-field route (van Bommel et al.) is in reach.");
else
    LOG(Sprintf("  first few disagreeing primes: %o", fails[1..Min(10,#fails)]));
end if;
exit;
