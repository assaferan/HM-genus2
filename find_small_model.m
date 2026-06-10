// find_small_model.m -- find shareable Weierstrass models for genus-2 curves.
//
// Key insight: Igusa-Clebsch invariants I2,I4,I6,I10 computed from a high-
// coefficient model have denominators ~10^63.  But I2 = 8*J2, I4 = 4*J2^2 -
// 96*J4, I6 = 8*J2^3 - 160*J2*J4 - 576*J6, I10 = 4096*J10 from the WPS-
// normalized CHIMP Igusa invariants gives IC with denominators <= 10^16.
// These small IC are WPS-equivalent (same isomorphism class).
// HyperellipticCurveFromIgusaClebsch on small IC gives a much smaller model.
//
// Run: magma find_small_model.m

AttachSpec("../CHIMP/CHIMP.spec");
import "Genus2Curve.m": Genus2CurveFromIgusa;

load "survivors_qsqrt5_gp11.m";   // defines k<s5>, survivors_qsqrt5_gp11
outfile := "small_models_qsqrt5_gp11.m";
System("rm -f " cat outfile);
PrintFile(outfile, "k<s5> := QuadraticField(5); small_models_qsqrt5_gp11 := [* *];");

function MaxCoeffSize(C)
    K := BaseRing(C); f, h := HyperellipticPolynomials(C);
    cs := [K!c : c in Coefficients(f) cat Coefficients(h)];
    cs := [c : c in cs | c ne 0];
    if #cs eq 0 then return 0; end if;
    return Max([Max([Abs(e) : e in Eltseq(c)]) : c in cs]);
end function;

function IntegralModelK(C)
    K := BaseRing(C); f, h := HyperellipticPolynomials(C);
    m := LCM([Integers()| Denominator(K!c) : c in Coefficients(f) cat Coefficients(h) cat [K!1]]);
    return HyperellipticCurve(Parent(f)!(m^2*f), Parent(f)!(m*h));
end function;

// Apply GL_2(K) transformation x -> (a*x+b)/(c*x+d) to a univariate sextic.
function ApplyGL2(f, a, b, c, d)
    K := BaseRing(Parent(f)); R := Parent(f); n := Degree(f);
    cs := [K| Coefficient(f, i) : i in [0..n]];
    result := R!0;
    for i in [0..n] do
        if cs[i+1] ne 0 then
            result +:= cs[i+1] * (R!(K!a)*R.1 + R!(K!b))^i
                               * (R!(K!c)*R.1 + R!(K!d))^(n-i);
        end if;
    end for;
    return result;
end function;

function ClearDenoms(f)
    K := BaseRing(Parent(f));
    m := LCM([Integers()| Denominator(K!c) : c in Coefficients(f) cat [K!1]]);
    return Parent(f)!(m * f);
end function;

// Compute IC = [I2,I4,I6,I10] from CHIMP's WPS-normalized J-invariants QI.
// Magma's conversion: I2=8*J2, I4=4*J2^2-96*J4, I6=8*J2^3-160*J2*J4-576*J6, I10=4096*J10.
// These IC have denominators <=max(den(QI)) instead of 10^63+ from C_init.
function SmallICFromQI(QI, K)
    J2 := K!QI[1]; J4 := K!QI[2]; J6 := K!QI[3]; J10 := K!QI[5];
    I2  := 8*J2;
    I4  := 4*J2^2 - 96*J4;
    I6  := 8*J2^3 - 160*J2*J4 - 576*J6;
    I10 := 4096*J10;
    return [K| I2, I4, I6, I10];
end function;

// Greedy GL_2(O_K) reduction: try a broad collection of small-entry transformations.
function ReduceByGL2(f0, K)
    phi := (1 + K.1)/2;   // golden ratio (unit in O_{Q(sqrt5)})
    C := HyperellipticCurve(f0);
    best_f := f0; best_sz := MaxCoeffSize(C);
    improved := true;
    // Keep iterating until no improvement.
    while improved do
        improved := false;
        // Integer shifts x -> x + b for b in [-8..8]
        for b in [-8..8] do
            for a in [1, -1] do
                try
                    f_new := ClearDenoms(ApplyGL2(best_f, K!a, K!b, K!0, K!1));
                    sz := MaxCoeffSize(HyperellipticCurve(f_new));
                    if sz lt best_sz then
                        best_sz := sz; best_f := f_new; improved := true;
                    end if;
                catch e; end try;
            end for;
        end for;
        // Inversions x -> 1/(x+n) and x -> (x+n)/x
        for n in [-5..5] do
            for tr in [[0,1,1,n], [1,n,1,0]] do
                try
                    f_new := ClearDenoms(ApplyGL2(best_f, K!tr[1], K!tr[2], K!tr[3], K!tr[4]));
                    sz := MaxCoeffSize(HyperellipticCurve(f_new));
                    if sz lt best_sz then
                        best_sz := sz; best_f := f_new; improved := true;
                    end if;
                catch e; end try;
            end for;
        end for;
        // O_K-transforms using golden ratio units.
        for u in [phi, phi^2, phi^(-1), phi^(-2)] do
            for v in [K|0, 1, -1, phi, -phi, phi^2, -phi^2] do
                for tr in [[u,v,K!0,K!1], [K!1,K!0,u,K!1], [u,v,K!1,K!0],
                           [u,K!0,K!0,K!1], [K!1,u,K!0,K!1]] do
                    try
                        f_new := ClearDenoms(ApplyGL2(best_f, tr[1], tr[2], tr[3], tr[4]));
                        sz := MaxCoeffSize(HyperellipticCurve(f_new));
                        if sz lt best_sz then
                            best_sz := sz; best_f := f_new; improved := true;
                        end if;
                    catch e; end try;
                end for;
            end for;
        end for;
    end while;
    return best_f, best_sz;
end function;

printf "Finding small models for %o survivor(s).\n\n", #survivors_qsqrt5_gp11;

for s in survivors_qsqrt5_gp11 do
    QI := s[2]; ht := s[3];
    printf "=== Ptors (height %o) ===\n", ht;
    printf "Compact Igusa invariants (for sharing):\n";
    printf "  J2=%o\n  J4=%o\n  J6=%o\n  J8=%o\n  J10=%o\n\n",
        QI[1], QI[2], QI[3], QI[4], QI[5];

    // Compute IC directly from QI (denominators <= 10^16, not 10^63+).
    IC := SmallICFromQI(QI, k);
    printf "Small IC: I2=%o\n  I4=%o\n  I6=%o\n  I10=%o\n\n",
        IC[1], IC[2], IC[3], IC[4];

    // Strategy 1: HyperellipticCurveFromIgusaClebsch on small IC.
    // This runs Mestre internally with HasRationalPoint on a small-coeff conic.
    best_C := 0; best_sz := -1; found := false;
    try
        C1 := HyperellipticCurveFromIgusaClebsch(IC);
        sz1 := MaxCoeffSize(C1);
        printf "Mestre (small IC): max coeff size %o\n", sz1;
        best_C := C1; best_sz := sz1; found := true;
    catch e;
        printf "Mestre (small IC) failed: %o\n", e`Object;
    end try;

    // Strategy 2: Build initial model from Genus2CurveFromIgusa and compare.
    try
        C_init := Genus2CurveFromIgusa(QI, k);
        sz_init := MaxCoeffSize(C_init);
        printf "Initial model (Genus2CurveFromIgusa): max coeff size %o\n", sz_init;
        if not found or sz_init lt best_sz then
            best_C := C_init; best_sz := sz_init; found := true;
        end if;
    catch e;
        printf "Genus2CurveFromIgusa failed: %o\n", e`Object;
    end try;

    if not found then
        printf "No model found for this survivor!\n\n";
        continue;
    end if;

    printf "Before GL_2: best max coeff size = %o\n\n", best_sz;

    // Strategy 3: GL_2(O_K) greedy reduction on integral model.
    CI := IntegralModelK(best_C);
    f0, _ := HyperellipticPolynomials(CI);
    f_red, sz_red := ReduceByGL2(f0, k);
    if sz_red lt MaxCoeffSize(best_C) then
        best_C := HyperellipticCurve(f_red);
        printf "GL_2 reduction improved to max coeff: %o\n", sz_red;
    else
        printf "GL_2 reduction: no improvement (max coeff %o)\n", MaxCoeffSize(best_C);
    end if;

    printf "\nBest model (max coeff %o):\n  %o\n\n", MaxCoeffSize(best_C), best_C;

    PrintFile(outfile, Sprintf("// ht=%o J10=%o", ht, QI[5]));
    PrintFile(outfile, Sprintf("Append(~small_models_qsqrt5_gp11, %o);", best_C));
end for;
printf "Written to %o.\n", outfile;
