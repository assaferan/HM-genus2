// The Horrocks--Mumford bundle: build it on P^4 and confirm its defining
// cohomology -- h^0(F) = 4 (four independent global sections) and h^0(F(-1)) = 0.
k := Rationals();                      // GF(p) with p prime is much faster;
P<[x]> := ProjectiveSpace(k, 4);
F := HorrocksMumfordBundle(P);
assert CohomologyDimension(F, 0,  0) eq 4;
assert CohomologyDimension(F, 0, -1) eq 0;
// The 4 degree-0 generators of the section module are the coordinates of the
// P^3 = P(H^0(F)) that parametrizes HM abelian surfaces.
M := GlobalSectionSubmodule(F);
secs := [M.i : i in [1..Rank(M)] | Grading(M)[i] eq 0];

// HMSurface: the (1,5)-polarized abelian surface A_x cut out by the section
// s = sum x_i * sec_i, for x in P^3. Returns the smooth degree-10 surface in P^4.
function HMSurface(F, x)
    M := GlobalSectionSubmodule(F);
    secs := [M.i : i in [1..Rank(M)] | Grading(M)[i] eq 0];
    assert #secs eq 4;
    s := &+[x[i]*secs[i] : i in [1..4]];      // the chosen global section
    Z := ZeroSubscheme(F, s);                 // its zero locus (the surface)
    I := Saturation(Ideal(Z));                // saturate to the true ideal
    return Scheme(Ambient(Z), I);
end function;

X := HMSurface(F, [1,0,3,7]);
assert Dimension(X) eq 2;
assert Degree(X) eq 10;
// The surface is cut out *set-theoretically* by 3 quintics, but the saturated
// ideal's minimal basis has more than 3 generators (extra higher-degree
// elements appear during saturation). So we filter for the degree-5 part
// instead of expecting #MinimalBasis eq 3 -- see x_to_tau in Inversion.m.
assert #[b : b in MinimalBasis(Ideal(X)) | Degree(b) eq 5] eq 3;
assert IsNonsingular(X);

// Complex uniformization: for A_tau = C^2 / (tau Z^2 + diag(1,5) Z^2),
// z --> ThetaPt(z, tau, N) is the (1,5) theta embedding A_tau --> P^4.
// The 5 coordinates are theta sums with characteristic offset (0, j/5), j = 0..4.

CC := ComplexField(200);
function ThetaPt(z, tau, N)            // z,tau over CC; returns P^4 point on A_tau
    ipi := Pi(CC)*CC.1; pt := [];
    for j in [0..4] do
        c := Matrix(CC,2,1,[0, j/5]); s := CC!0;
        // Truncated lattice sum over |n1|,|n2| <= N; converges as Im(tau) > 0.
        for n1 in [-N..N], n2 in [-N..N] do
            v := Matrix(CC,2,1,[n1,n2]) + c;
            // exp(pi i v^t tau v + 2 pi i v^t z)
            s +:= Exp(ipi*(Transpose(v)*tau*v)[1,1] + 2*ipi*(Transpose(v)*z)[1,1]);
        end for;
        Append(~pt, s);
    end for;
    return pt;
end function;