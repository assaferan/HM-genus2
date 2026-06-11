# Project handoff — HM-genus2 (and the new (1,11) extension)

Status snapshot for picking up on another machine. Repo: `assaferan/HM-genus2` (GitHub,
public, MIT). Branch `main`, all work committed and pushed. This file is the orientation;
`CLAUDE.md` and `README.md` are the canonical per-file docs.

## What this project is

Two parallel pipelines that reconstruct a genus-2 curve from a point parametrizing an
abelian surface, so that `Jac(C)` carries a **rational cyclic isogeny**:

1. **(1,5) — Horrocks–Mumford (DONE, production).**
   `x ∈ P³ ──HM bundle──► A_x ⊂ P⁴ (1,5)-polarized ──invert theta──► τ ──► genus-2 C`
   with a rational **5-isogeny**. Fully working, including a number-field version and a
   search for small-conductor curves over `ℚ(√2)`.

2. **(1,11) — Gross–Popescu (NEW, in progress).**
   `v (Klein-cubic moduli) ──► A_v ⊂ P¹⁰ (1,11)-polarized ──invert theta──► τ ──► genus-2 C`
   with a rational **11-isogeny**. The hard parts are verified and the inversion converges
   (see below); what remains is wrapping `τ → curve` and building a search.

## Environment / dependencies

- **Magma** (the whole thing is Magma). Invoked as `magma file.m`.
- **CHIMP** must be checked out as a **sibling directory** `../CHIMP` (so
  `AttachSpec("../CHIMP/CHIMP.spec")` resolves). Provides `ComplexFieldExtra`,
  `ReduceSmallPeriodMatrix`, `GeometricEndomorphismRepresentationCC`,
  `HyperellipticCurveFromIgusaInvariants`, etc.
- `*.sig` are Magma attach caches (ignored). Scratch scripts are `_*.m` (gitignored);
  `survivors.m`, `search_results.txt` are run outputs (gitignored).

## Files (committed)

**(1,5) pipeline — works:**
- `HMsurface.m` — HM bundle; `HMSurface(F,x)`; `ThetaPt`.
- `Inversion.m` — reference inverse (finite-diff Gauss–Newton) + `verify_tau`.
- `InversionFast.m` — **production** inverse `x_to_tau_fast`: analytic Jacobian
  (`ThetaJet`), Levenberg–Marquardt (`FindTauFast`), adaptive truncation, two-phase
  precision, embedding round-robin for number fields. **Read this — the (1,11) inversion
  is a direct port of it.**
- `Genus2Curve.m` — `Genus2Curve(x : K:=...)`: `x → τ → curve/K`. Also
  `Genus2CurveFromTau`, `Genus2CurveFromIgusa(QI,K)` (descend from invariants, no period
  recomputation), `IgusaNumericFromTau`, `IgusaInvariantsInK`, and
  `GeometricEndomorphismDimension(τ)` (=1 ⟺ `End(J)=ℤ`).
- `Reduction.m` — `PotentialGoodReductionAt5`, `GoodReductionAt5`, `ConductorExponentAt`,
  `ReductionReport`, `MinimalTwist` (twist to good reduction at every prime of potential
  good reduction → smallest conductor), `FrobeniusPolynomial`, `IntegralModelOf`.
- `Heights.m` — `BoundedHeightPoints(K, M : ExcludeRational)`; `IsDefinedOverQ`.
- `search_quadratic.m` — **screening** (fast): enumerate non-rational `ℚ(√2)` points by
  height; keep genuine-over-`K` (not a base change), potential-good-at-5, simple
  (`End=ℤ`); stream survivors + their Igusa invariants to `survivors.m`.
- `analyze_survivors.m` — **analysis** (slow, on-demand): from `survivors.m`, build curve,
  `MinimalTwist`, conductor exponents, factored Frobenius-at-5. (Split out because the
  descent + ~30 conductor computations per survivor are expensive — never put them in the
  per-point screening loop.)
- `test.m`, `test_quadratic.m`.

**(1,11) pipeline — new:**
- `GP11.m` — the explicit Gross–Popescu construction (see below). Has `ThetaPt11`,
  `HeisenbergGenerators11`, `KleinCubicMatrix`, `IntertwiningMatrix`,
  `QuadricsFromVector/Pencil`, `SurfaceQuadrics`, `PencilFromTorsion`, `GP11SelfTest`,
  and an EARLY finite-difference inversion (`ResidualAtZ11`, `FindTau11`, `InvertGP11`) —
  **superseded** by `InversionGP11.m`; keep using the latter.
- `InversionGP11.m` — **the working (1,11) inversion** (analytic-Jacobian port of
  `InversionFast.m`): `ThetaJet11`, `QuadricsCC`, `ResJac11`, `FindTauFast11`,
  `InvertGP11Fast`.

## The (1,11) Gross–Popescu construction (verified)

Refs: **GP2** = Gross–Popescu, *The moduli space of (1,11)-polarized abelian surfaces is
unirational*, arXiv **math/9902017** (the source of the formulas below); **GP1** =
*Equations of (1,d)-polarized abelian surfaces*, Math. Ann. 310 (1998); Adler, Amer. J.
Math. 100 (1978). Moduli `A₁₁^lev` is birational to the **Klein cubic**
`K = {x₀²x₁+x₁²x₂+x₂²x₃+x₃²x₄+x₄²x₀=0} ⊂ P⁴`.

- **Heisenberg `H₁₁` on `P¹⁰`** (coords `x₀..x₁₀`): `σ(x_i)=x_{i-1}`, `τ(x_i)=ξ⁻ⁱx_i`
  (`ξ=e^{2πi/11}`), `ι(x_i)=x_{-i}`.
- **Quadrics**: `S²(V)` (66-dim) = 6 copies of an 11-dim `H₁₁`-rep, organized by
  `(R5)_{ij}=x_{j+i}x_{j-i}` (`0≤i≤5, 0≤j≤10`, mod 11). A vector `v=(v₀..v₅) ∈ V₊≅ℂ⁶`
  gives 11 quadrics `Q^v_j = Σ_i v_i x_{j+i}x_{j-i}`. **A (1,11) surface = the 22 quadrics
  of a pencil `⟨v,w⟩ ⊂ V₊`** (a point of `Gr(2,6)`; `im(Θ₁₁) ⊂ Gr(2,6)` is birational to
  the Klein cubic, cut out by 5 linear Plücker relations — see `PluckerRelationsX`).
- **Two 6×6 skew matrices** (in `GP11.m`): `M` (linear in `x₀..x₄`), `Pf(M)=`Klein cubic;
  and the **intertwining `S`** (quadratic in `x₁..x₅`), `Pf(S)=f₆` the sextic `D₂⊂P⁴` of
  odd 2-torsion points.
- **The explicit route to a surface (use this):** `P ∈ D₂  ↦  ker S(P) = pencil ⟨v,w⟩  ↦
  22 quadrics` (`SurfaceQuadrics`/`PencilFromTorsion`). The Klein-cubic↔`Gr(2,6)` map is
  the *indirect* Fano–Iskovskih one — don't bother; `D₂` (or `X⊂Gr(2,6)`) is direct.

### Verified facts (don't re-litigate)
- `Pf(M) = Klein cubic`, `Pf(S) =` sextic, pencil → 22 deg-2 quadrics. (`GP11SelfTest`.)
- **Convention check passed:** for a random τ, the (1,11) theta points lie on **exactly a
  2-dimensional pencil** of `R5`-quadrics (the `M^H M` of the `R5`-evaluation matrix has
  two machine-zero eigenvalues, cleanly separated). So `ThetaPt11` ↔ GP quadrics are
  consistent — the construction is correct end-to-end.
- **The inversion converges:** from a real `D₂` point, `InvertGP11Fast` reached
  `|R| = 8.85e-59` at precision 60 (basin hit rate ≈ 1/13 ≈ 8%, like (1,5)).

### Critical gotchas for (1,11) (these cost time to discover)
1. **Use a REAL root** of the `D₂`-defining polynomial. A generic `D₂` point: fix
   `x₁..x₄` rational (e.g. `1,2,3,4`), set `x₅ = ` a **real** root of `Pf(S)(1,2,3,4,t)`
   (for `1,2,3,4` this is the irreducible quartic `3t⁴+16t³+124t²+174t−1190`, 2 real
   roots ≈ `−4.184, 2.169`). A **complex** embedding gives a non-real surface with no τ
   in the standard domain → inversion never converges.
2. The `D₂` point must have `Rank(S)=4` (kernel dim 2; not on `D₁`). Small-integer points
   like `[-3,-3,0,0,0]` are **degenerate** (too many zero coords) — avoid.
3. The inversion residual must use the **raw** theta point with an **analytic** Jacobian
   that dehomogenizes at the largest coordinate (chart fixed *within* an evaluation). A
   *finite-difference* Jacobian with max-coordinate normalization is BROKEN (the argmax
   switches between base and perturbed points → garbage Jacobian → stuck at `|R|≈√3`).
   `InversionGP11.m` does it right; the early `FindTau11` in `GP11.m` does not.

## What remains for (1,11)

1. **Wrap `τ → curve` — CONFIRMED working.** `A_τ` (1,11-polarized) is 11-isogenous to the
   principally polarized `[τ|I]`, so the **same τ** from `InvertGP11Fast` feeds the reused
   `IgusaNumericFromTau` / `Genus2CurveFromTau` (in `Genus2Curve.m`). The full-chain test
   ran successfully: from the real `D₂` point above, `|R|=9.95e-79` and the genus-2 curve's
   absolute Igusa invariants came out real (curve over a number field):
   `[1, -0.158469023628550, -0.188271294021765, -0.0533459313678878, -2.0704e-7]`.
   So the entire chain `v → A ⊂ P¹⁰ → τ → C` works. What's left here is packaging it as a
   tidy `GP11Curve(...)` driver and recognizing the invariants in the field of moduli.
2. **11-isogeny test** — analogue of `test.m`: at good primes `p ≠ 11`, the mod-11
   reduction of the `L`-polynomial of `C/𝔽_p` has a root in `𝔽₁₁`.
3. **Search infrastructure** — analogue of `Heights.m` + `search_quadratic.m`: enumerate
   rational points on the Klein cubic (or on `D₂`), build surfaces, invert, reconstruct,
   filter by good-reduction-at-11 + small conductor + `End=ℤ` + not-a-base-change. NB the
   easy rational `D₂` points are degenerate; finding `D₂` points whose **field of moduli is
   small** (so the curve is over ℚ or a small field) is itself a search.
4. Reuse `Reduction.m`/`Genus2Curve.m` with "5 → 11" throughout (good reduction at 11,
   `L_p mod 11`, minimal twist, conductor — all polarization-degree-agnostic).

### Reproduce the (1,11) end-to-end run
```magma
AttachSpec("../CHIMP/CHIMP.spec");
import "GP11.m": IntertwiningMatrix, PencilFromTorsion;
import "InversionGP11.m": InvertGP11Fast;
import "Genus2Curve.m": IgusaNumericFromTau;
Qt<t> := PolynomialRing(Rationals());
f := Pfaffian(IntertwiningMatrix([Qt| 1,2,3,4,t]));           // = 3t^4+16t^3+124t^2+174t-1190
g := [h[1] : h in Factorization(f) | Degree(h[1]) ge 1][1];   // the irreducible quartic
K<a> := NumberField(g);
v, w := PencilFromTorsion([K| 1,2,3,4,a]);
CC := ComplexFieldExtra(80);
rt := CC ! [r[1] : r in Roots(g, RealField(80))][1];          // a REAL root (~ -4.184)
emb := hom< K -> CC | rt >;
tau, nr := InvertGP11Fast([emb(c):c in v], [emb(c):c in w], CC : trials := 80);  // |R| -> ~0
J := IgusaNumericFromTau(tau);                                // genus-2 curve invariants
```

## (1,5) search — current state

- **27 survivors** found over `ℚ(√2)` (BoxM=2, MaxPts=30, Trials=120); stored in
  `survivors.m`. Extended search (`search_quadratic_ext.m`, MaxPts=300, Trials=300)
  is available to run for more.
- All survivors have good reduction at 5 (both primes above 5) and simple Jacobian.
- Full Frobenius data in `data_qsqrt2.m` / `data_qsqrt2.txt`.
- **Known issue:** `TryReduceModP` in `data_qsqrt2.m` fails at both primes above 3 and 5
  for all 27 survivors. This is a **mathematical obstruction** (not a bug): the Mestre
  reconstruction gives a model where all P₅-valuations of the polynomial coefficients
  are negative (because v_{P5}(J₁₀)/10 requires a scaling u with v_{P5}(u) = 1/2,
  which is not achievable over K). Working at these primes requires passing to the
  ramified extension K(5^{1/2}).

## (1,11) search — current state (NEW, COMPLETE)

- **3 K-isomorphism classes** found over `ℚ(√5)` (BoxM=2, MaxPts=200): Curve A, B, B'.
  Stored in `survivors_qsqrt5_gp11.m`.
- **Curve A**: smallest conductor (odd part = norm 66179). Fully verified:
  - Bad prime N(p)=66179, conductor exponent 1 ✓
  - Both primes above 11: conductor exponent 0 (good reduction) ✓
  - Frobenius L-polynomials verified against `data_qsqrt5_gp11.txt` ✓
  - Verification script: `verify_curveA.m` (~8 s)
- **Curve B / B'**: Galois conjugate pair, conductor exponent 1 at N(p)=602479.
- Full Frobenius data: `data_qsqrt5_gp11.txt`.  Results summary: `results_gp11.md`.
- Extended search: `search_qsqrt5_gp11_ext.m`.
- **Weierstrass model note:** coefficients have 170–350 digit numerators/denominators;
  this is intrinsic (see `results_gp11.md` for explanation and `reduce_curveA*.m`
  for the reduction attempts). For computation use `TryReduceModP`; for a paper
  quote the Igusa-Siegel invariants and conductor.

## Git
All committed and pushed to `main`. Recent commits:
- `65177eb` — extended searches, `verify_curveA.m`, improved `TryReduceModP`.
- Earlier: `InversionGP11.m`, GP11 convention verification, search infrastructure.
`git pull` gets everything; remember to also have `../CHIMP` as a sibling directory.
