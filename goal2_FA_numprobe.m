SetColumns(0);
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
C := HyperellipticCurve(-10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9);
prec := 100;
AJ := AnalyticJacobian(C : Precision := prec);
BPM := BigPeriodMatrix(AJ);      // 2 x 4, columns = lattice basis
printf "period matrix %o x %o at precision %o\n", Nrows(BPM), Ncols(BPM), prec;
// a 7-torsion point: z = (1/7) * (column of the period lattice) mod Lambda
// take lattice vector v = BPM[.,1..4] combos; simplest nonzero 7-torsion:
CC := BaseRing(BPM);
lat := [ Matrix(CC,2,1,[BPM[1,j],BPM[2,j]]) : j in [1..4] ];
// probe: what maps a point on the analytic Jacobian to algebraic coordinates?
printf "AnalyticJacobian type: %o\n", Type(AJ);
// list relevant intrinsics
z := (1/7)*lat[1];              // a 2x1 complex vector, a 7-torsion point
printf "sample 7-torsion z computed (first coord ~ %o)\n", z[1,1];
try
    pt := FromAnalyticJacobian(z, AJ);
    printf "FromAnalyticJacobian OK: image is a sequence of %o points on C\n", #pt;
    printf "  first point coords: %o\n", pt[1];
catch e
    printf "FromAnalyticJacobian failed/na: %o\n", e`Object;
end try;
exit;
