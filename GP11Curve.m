// Driver for the (1,11) / Gross-Popescu pipeline:
//   D_2 point Ptors (over K, with embedding emb: K -> CC)
//   -> pencil <v,w> -> tau -> genus-2 curve C with a rational 11-isogeny.
//
// Usage:
//   Qt<t> := PolynomialRing(Rationals());
//   f := Pfaffian(IntertwiningMatrix([Qt| 1,2,3,4,t]));
//   g := [h[1] : h in Factorization(f) | Degree(h[1]) ge 1][1];
//   K<a> := NumberField(g);
//   CC := ComplexFieldExtra(120);
//   rt := CC ! [r[1] : r in Roots(g, RealField(120))][1];  // pick a REAL root
//   emb := hom< K -> CC | rt >;
//   C, QI, tau := GP11Curve([K| 1,2,3,4,a], K, emb, CC);
AttachSpec("../CHIMP/CHIMP.spec");
import "GP11.m": PencilFromTorsion;
import "InversionGP11.m": InvertGP11Fast;
import "Genus2Curve.m": Genus2CurveFromTau;

// Main driver for the (1,11) pipeline.
// Ptors: D_2 point [x1..x5] over K; emb: K -> CC; CC: ComplexFieldExtra.
// Returns (curve C over K, absolute Igusa invariants in K, period matrix tau).
function GP11Curve(Ptors, K, emb, CC : Trials := 100, Verbose := false)
    v, w := PencilFromTorsion(Ptors);
    vCC := [emb(c) : c in v];
    wCC := [emb(c) : c in w];
    tau, nr := InvertGP11Fast(vCC, wCC, CC : trials := Trials, verbose := Verbose);
    if Verbose then printf "InvertGP11Fast: |R| = %o\n", RealField(8)!nr; end if;
    C, QI := Genus2CurveFromTau(tau, K, emb);
    return C, QI, tau;
end function;
