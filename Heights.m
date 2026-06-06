// Enumerate K-rational points of P^3 (K a real quadratic field) of bounded height.
//
// Each point has a primitive O_K representative (coordinates with trivial gcd; K has
// class number 1 for the fields of interest). We enumerate representatives with
// coordinates a + b*w (w = K.1) in the box |a|,|b| <= M, deduplicate each projective
// point by its canonical affine form (divide by the LAST nonzero coordinate, which
// removes all scaling -- including unit -- ambiguity), and sort by the multiplicative
// Weil height H = prod_v max_i |x_i|_v (v the archimedean places; the non-archimedean
// part is 1 for a primitive representative). Optionally keep only H <= HeightBound.
//
// Note: to capture *every* point of height <= B, take M >= B (a height-B point's
// primitive coordinates are bounded by ~B in each embedding); for the small B that the
// per-point period computation allows, the box is small.
function BoundedHeightPoints(K, M : HeightBound := 0)
    Ok := Integers(K);
    w := K.1;
    box := [ K | a + b*w : a in [-M..M], b in [-M..M] ];
    seen := { };
    res := [];
    for tup in CartesianPower(box, 4) do
        x := [ tup[1], tup[2], tup[3], tup[4] ];
        nz := [ i : i in [1..4] | x[i] ne 0 ];
        if #nz eq 0 then continue; end if;
        j := nz[#nz];                              // canonical affine chart: last nonzero = 1
        key := < x[i]/x[j] : i in [1..4] >;
        if key in seen then continue; end if;
        Include(~seen, key);
        g := GCD([ Ok | x[i] : i in nz ]);         // make primitive (class number 1)
        xp := [ K | x[i]/g : i in [1..4] ];
        cs := [ Conjugates(c) : c in xp ];         // cs[i] = [ sigma_1(x_i), sigma_2(x_i) ]
        H := &*[ Max([ Abs(cs[i][t]) : i in [1..4] ]) : t in [1..#cs[1]] ];
        Append(~res, < xp, H >);
    end for;
    flt := [ r : r in res | HeightBound eq 0 or r[2] le HeightBound ];
    Sort(~flt, func< a, b | a[2] lt b[2] select -1 else (a[2] gt b[2] select 1 else 0) >);
    return [ r[1] : r in flt ], [ RealField(6)!r[2] : r in flt ];
end function;
