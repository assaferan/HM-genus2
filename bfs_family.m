// bfs_family.m -- Bruin-Flynn-Shnidman explicit family of genus-2 curves with
// RM by sqrt(3) (arXiv:2102.04319, Section 3, "sqrt3 over the ground field").
//
//   C_{a,b,c}: y^2 = G1(x)^2 + lam1*H1(x)^3   (= q4^2 * f_{a,b,c}, iso to y^2=f).
//
// Validate the transcription against the paper's sample C_{1,2,-1}, whose
// Jacobian (Thm 1.4) is y^2 = 8x^5 - 3x^4 - 2x^3 - 7x^2 + 4x + 20.
//
// Run: magma bfs_family.m
SetColumns(0);

// Q(X,Y,Z) = X^2 + X(Y-Z) + (Y-Z)^2 - 3XZ   (eq (7))
function Qform(X,Y,Z) return X^2 + X*(Y-Z) + (Y-Z)^2 - 3*X*Z; end function;

// t3 = Q(a, zeta b, c) * Q(a, zeta^2 b, c), rational in a,b,c.  Compute via zeta.
K3<z3> := CyclotomicField(3);
function t3val(a,b,c)
    v := Qform(a, z3*b, c) * Qform(a, z3^2*b, c);   // lands in the base field
    return v;
end function;

// The family sextic for (a,b,c) over a field containing the base.  Returns the
// degree-<=6 polynomial G1^2 + lam1*H1^3 in Px.
function BFSsextic(a, b, c, Px)
    x := Px.1;
    q4 := Qform(c,b,a);
    t3 := t3val(a,b,c);
    // coerce t3 (base-field value stored in K3-compositum) back to the coeff ring
    if Type(t3) eq FldCycElt then t3 := Rationals()!t3; end if;
    H1  := b*(a-c)*q4*x^2 - q4*x - 1;
    lam1 := 4*a*c*t3;
    quart := a^4 - 4*a^3*b - 8*a^3*c + 3*a^2*b^2 + 18*a^2*b*c + 18*a^2*c^2
             + 2*a*b^3 - 9*a*b^2*c - 12*a*b*c^2 - 8*a*c^3 - 2*b^4 - b^3*c
             + 2*b*c^3 + c^4;
    cub := 2*a^3 - 3*a^2*b - 9*a^2*c + 9*a*b*c + 6*a*c^2 + b^3 - c^3;
    con := a^2 - 2*a*b - 4*a*c + b^2 + 7*b*c + c^2;
    G1 := b*(a-c)^2*(a^2-4*a*c-b^2+c^2)*q4^2*x^3
          - (a-c)*q4*quart*x^2
          - q4*cub*x
          - a*con;
    return G1^2 + lam1*H1^3;
end function;

QQ := Rationals(); Px<x> := PolynomialRing(QQ);

// --- validation ---
fam := BFSsextic(1, 2, -1, Px);
print "family sextic C_{1,2,-1} (= q4^2 f):", fam;
Cfam := HyperellipticCurve(fam);
Csamp := HyperellipticCurve(8*x^5 - 3*x^4 - 2*x^3 - 7*x^2 + 4*x + 20);

igF := IgusaClebschInvariants(Cfam);
igS := IgusaClebschInvariants(Csamp);
print "Igusa-Clebsch (family) :", igF;
print "Igusa-Clebsch (sample) :", igS;

// projective comparison (weighted): normalize by first nonzero ratio
function ProjEq(u, v)
    // compare weighted-projective tuples up to a common scalar^weight; just test
    // all cross ratios I_i^{w_j} vs I_j^{w_i} are equal via pairwise scaling.
    // Simple robust test: exists t with u[i] = t^{w_i} v[i]; use weights (2,4,6,10).
    w := [1,2,3,5];   // Igusa-Clebsch weights /2
    // find scale from first nonzero
    for i in [1..4] do
        if v[i] ne 0 and u[i] ne 0 then
            // t^{w_i} = u[i]/v[i]
            for j in [1..4] do
                if u[i]^w[j]*v[j]^w[i] ne v[i]^w[j]*u[j]^w[i] then return false; end if;
            end for;
            return true;
        end if;
    end for;
    return false;
end function;

print "ISOMORPHIC (projective Igusa-Clebsch match):", ProjEq(igF, igS);
print "IsIsomorphic:", IsIsomorphic(Cfam, Csamp);
exit;
