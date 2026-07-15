// goal2_mod5_idx5_cert2.m -- resume the GRH trace certificate for idx 5 from
// where goal2_mod5_idx5.m was cut off by the 55-min timeout (l ~ 5003).
// Verifies good primes in [5000, 12000], appending to the tally.
SetColumns(0);
d := 2; p := 5;
F<w> := QuadraticField(d); OF := Integers(F);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
fc := [5,-10,10,0,5,2];
fpol := &+[ fc[i]*x^(i-1) : i in [1..#fc] ];
C := HyperellipticCurve(fpol);
Fp := GF(p); badp := {2,5};

LOGF := "goal2_mod5_idx5_cert2_out.txt";
PrintFile(LOGF, "# cert2: resume GRH check for idx 5, primes [5000,12000]" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

P2 := Factorization(2*OF)[1][1]; P5 := Factorization(5*OF)[1][1];
nlev := P2^3 * P5^2; assert Norm(nlev) eq 5000;
M := HilbertCuspForms(F, nlev, [2,2]);
Dc := NewformDecomposition(NewSubspace(M));
ef := Eigenform(Dc[6]); E := HeckeEigenvalueField(Dc[6]); OE := Integers(E);
found := false; Fq := Fp; red := func< a | Fp!(Integers()!a) >;
for pl in Factorization(5*OE) do
    Fq0, r0 := ResidueClassField(pl[1]);
    Fq := Fq0; red := func< a | r0(OE!a) >; found := true; break;
end for;
assert found;
LOG(Sprintf("orbit-6 eigenform ready, residue F_%o. Verifying [5000,12000] ...", #Fq));

LO := 5000; HI := 12000;
gsplit := 0; ginert := 0; fails := []; nextlog := LO + 1000;
for l in PrimesInInterval(LO, HI) do
    if l in badp then continue; end if;
    aA := Fp ! (l + 1 - #ChangeRing(C, GF(l)));
    if IsSquare(GF(l)!d) then
        fp := Factorization(l*OF);
        s := red(OE!HeckeEigenvalue(ef, fp[1][1])) + red(OE!HeckeEigenvalue(ef, fp[2][1]));
        ok := (s eq Fq!aA); gsplit +:= 1;
    else
        ok := (aA eq Fp!0); ginert +:= 1;
    end if;
    if not ok then Append(~fails, l); end if;
    if l gt nextlog then
        LOG(Sprintf("  ... up to %o: %o split + %o inert checked, %o fails", l, gsplit, ginert, #fails));
        nextlog +:= 1000;
    end if;
end for;
LOG(Sprintf("\n[5000,12000]: %o good primes (%o split, %o inert), %o disagreements",
    gsplit+ginert, gsplit, ginert, #fails));
if #fails ne 0 then LOG(Sprintf("  fails: %o", fails[1..Min(10,#fails)])); end if;
exit;
