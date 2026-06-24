// prescan_rational_gp11.m -- Cheap (no-inversion) hunt for RATIONAL non-degenerate D_2
// points, which give (1,11) abelian surfaces / genus-2 curves over Q.
// D_2 point = [a1,a2,a3,a4,x5] with x5 a RATIONAL root of g(t)=Pf(S)(a1,a2,a3,a4,t),
// and Rank(S)=4 (kernel dim 2; non-degenerate, not on D_1).
//
// Run: magma prescan_rational_gp11.m

SetColumns(0);
import "GP11.m": IntertwiningMatrix;
Qt<t> := PolynomialRing(Rationals());

BoxM := 8;     // |a_i| <= BoxM
found := [];

function NonDegenerate(P)   // P a 5-vector over Q; need rank S = 4
    S := IntertwiningMatrix(P);
    return Rank(S) eq 4;
end function;

count := 0;
for a1 in [-BoxM..BoxM], a2 in [-BoxM..BoxM], a3 in [-BoxM..BoxM], a4 in [-BoxM..BoxM] do
    // primitive-ish dedup: require gcd 1 and first nonzero > 0
    cs := [a1,a2,a3,a4];
    if cs eq [0,0,0,0] then continue; end if;
    g0 := Gcd([Integers()| c : c in cs]);
    if g0 ne 1 then continue; end if;
    fnz := [c : c in cs | c ne 0][1];
    if fnz lt 0 then continue; end if;
    g := Pfaffian(IntertwiningMatrix([Qt| a1,a2,a3,a4,t]));
    if g eq 0 then continue; end if;
    for r in Roots(g) do
        t0 := r[1];
        P := [Rationals()| a1,a2,a3,a4,t0];
        // avoid very degenerate (too many zero coords) and require rank 4
        if #[c : c in P | c eq 0] le 2 and NonDegenerate(P) then
            count +:= 1;
            Append(~found, P);
        end if;
    end for;
end for;

printf "Rational non-degenerate D_2 points with |a_i|<=%o: %o\n\n", BoxM, #found;
// sort by a naive height and show the smallest
ht := func<P | Maximum([Abs(Numerator(c)) : c in P] cat [Abs(Denominator(c)) : c in P])>;
Sort(~found, func<x,y | ht(x)-ht(y)>);
for i in [1..Min(25,#found)] do
    printf "  %o   (rank S=4, height %o)\n", found[i], ht(found[i]);
end for;

// write them for the inversion stage
PrintFile("rational_d2_points.txt", "// rational non-degenerate D_2 points [a1..a4,x5]" : Overwrite := true);
PrintFile("rational_d2_points.txt", Sprintf("d2pts := %o;", found));
printf "\nWrote rational_d2_points.txt (%o points).\n", #found;
