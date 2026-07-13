// mod7_twist_check.m -- are the 5 mod-7 examples (idx 33,53,54,73,74) the same
// residual rep up to quadratic twist?  a_l^2 mod 7 is quadratic-twist-invariant;
// if it agrees across all five they are twists of one form (the degree-18 one).
SetColumns(0);
load "goal1_data.m";
Px<x> := PolynomialRing(Rationals()); F7:=GF(7);
idxs := [33,53,54,73,74];
curves := [];
for i in idxs do
    e := [z:z in entries|z[1] eq i][1]; fc:=e[5]; hc:=e[6];
    Append(~curves, HyperellipticCurve(&+[fc[j]*x^(j-1):j in [1..#fc]], &+[hc[j]*x^(j-1):j in [1..#hc]]));
end for;
printf "  l  | a_l^2 mod 7  for idx 33,53,54,73,74 | same?\n";
allsame := true;
for l in PrimesInInterval(11,120) do
    if l in {2,3,5,7} then continue; end if;
    sq := [ Integers()!((F7!(l+1-#ChangeRing(C,GF(l))))^2) : C in curves ];
    if #Seqset(sq) ne 1 then allsame := false; end if;
    printf " %3o | %o | %o\n", l, sq, #Seqset(sq) eq 1 select "yes" else "NO";
end for;
printf "\nAll five are quadratic twists of one residual rep: %o\n", allsame;
exit;
