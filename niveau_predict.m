SetColumns(0);
load "goal1_data.m";
Px<x> := PolynomialRing(Rationals());
F5 := GF(5); P5<t> := PolynomialRing(F5);
// multiplicity partition of the 6 roots of 4f+h^2 mod 5 (roots at infinity when deg<6)
function Niveau(fc, hc)
    f := &+[fc[i]*x^(i-1):i in [1..#fc]]; h := &+[hc[i]*x^(i-1):i in [1..#hc]];
    G := P5!(4*f + h^2);
    if G eq 0 then return "N2", [6]; end if;   // non-reduced / fully degenerate
    part := [ a[2] : a in Factorization(G) ];   // finite root multiplicities
    if Degree(G) lt 6 then Append(~part, 6 - Degree(G)); end if;  // roots at infinity
    Sort(~part); Reverse(~part);
    // rule: contains a 5 -> niveau-1 ; else -> niveau-2
    niv := (5 in part) select "N1" else "N2";
    return niv, part;
end function;
PrintFile("niveau_predict.txt","# idx : predicted-niveau : root-multiplicity-partition at 5" : Overwrite);
for ent in entries do
    idx := ent[1]; if ent[3] ne 5 then continue; end if;
    niv, part := Niveau(ent[5], ent[6]);
    PrintFile("niveau_predict.txt", Sprintf("%o %o %o", idx, niv, part));
    printf "idx %2o -> %o  %o\n", idx, niv, part;
end for;
exit;
