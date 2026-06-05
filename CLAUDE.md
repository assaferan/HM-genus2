# CLAUDE.md

Guidance for working with this repository.

## Overview

This is a **Magma** (computational algebra) project that reconstructs a genus-2
curve from a point parametrizing an abelian surface via the **Horrocks–Mumford
(HM) bundle**. The end-to-end pipeline is:

```
x ∈ P³  ──HMSurface──►   abelian surface A_x ⊂ P⁴   (smooth, degree 10, (1,5)-polarized)
        ──x_to_tau──►    period matrix τ ∈ H₂      (2×2 symmetric, Im τ > 0)
        ──Genus2Curve──► genus-2 curve C with Jac(C) ≅ A_τ
```

The forward theta embedding `A_τ → P⁴` is transcendental and has no closed-form
inverse, so `x_to_tau` recovers `τ` **numerically**: it evaluates the surface's
defining quintics at theta images of random points and solves the resulting
least-squares problem for the three entries of `τ`.

## Mathematical background

- **Horrocks–Mumford bundle** `F`: the unique (up to twist) indecomposable
  rank-2 vector bundle on `P⁴` with `h⁰(F) = 4`. Its sections vanish on abelian
  surfaces. A choice of section `s = Σ xᵢ·secᵢ` (with `x ∈ P³`, since
  `h⁰(F) = 4`) cuts out a `(1,5)`-polarized abelian surface `A_x ⊂ P⁴` of
  degree 10.
- **(1,5) theta embedding**: for `τ` in the genus-2 Siegel upper half-space,
  `A_τ = C² / (τ Z² + diag(1,5) Z²)`, and `z ↦ ThetaPt(z, τ)` embeds `A_τ` into
  `P⁴` using the 5 theta functions with characteristics `j/5`, `j = 0..4`.
- **Inversion principle**: a point `z` lies on `A_τ`, so its image `Θ(z)` must
  satisfy every defining equation `q` of the surface `A_x`. Treating
  `q(Θ(z))` as a residual that depends on `τ = [[a,b],[b,c]]`, we solve for
  `(a,b,c)` by nonlinear least squares over several random sample points `z`.

## Files

| File              | Role |
|-------------------|------|
| `HMsurface.m`     | Builds the HM bundle and the abelian surface `A_x = HMSurface(F, x)`; defines the reference `ThetaPt` embedding. |
| `Inversion.m`     | Reference inverse `x_to_tau`: finite-difference Gauss–Newton (`FindTau`) with random restarts and verification (`verify_tau`). |
| `InversionFast.m` | Optimized inverse `x_to_tau_fast`: analytic Jacobian (`ThetaJet`), Levenberg–Marquardt (`FindTauFast`), adaptive lattice truncation. **This is the one used in production.** |
| `Genus2Curve.m`   | Top-level driver: `x → τ → curve` via CHIMP's `ReconstructCurve`. |
| `test.m`          | Smoke test: from `x = [1,2,3,4]`, check the Jacobian has 5-torsion. |
| `*.sig`           | Magma-generated attach caches — **not source**, safe to ignore/regenerate. |

## Key functions

- `HMSurface(F, x)` — global section `s = Σ xᵢ secᵢ`, returns `Scheme` of the
  saturated zero ideal of `s`. Surface is dimension 2, degree 10, nonsingular.
- `ThetaPt(z, a, b, c, N, CC)` — truncated `(1,5)` theta point of `z`, lattice
  sum over `|n₁|,|n₂| ≤ N`. (`HMsurface.m` has a variant taking `tau` as a
  matrix.)
- `ThetaJet(...)` (in `InversionFast.m`) — returns theta values **and** their
  three `τ`-derivatives in a single lattice pass; the analytic Jacobian.
- `FindTau` / `FindTauFast` — least-squares solve for `(a,b,c)`. Fast variant
  uses Levenberg–Marquardt damping, backtracking, positive-definiteness guards,
  and ramps the truncation `N` up only as the residual shrinks.
- `x_to_tau` / `x_to_tau_fast` — random restarts, keeps best positive-definite
  solution, verifies on fresh points and the theta-null, returns the symmetric
  `τ`.
- `Genus2Curve(x : Prec := 200)` — composes the above with CHIMP's
  `ReconstructCurve`.

## Dependencies

- **Magma** (the HM bundle, theta sums, scheme/ideal machinery are all Magma).
- **CHIMP** — attached in `Genus2Curve.m` via
  `AttachSpec("../CHIMP/CHIMP.spec")`. Provides `ComplexFieldExtra`,
  `RationalsExtra`, and `ReconstructCurve`. **CHIMP must be checked out as a
  sibling directory** (`../CHIMP`).

## Running

```magma
// Build a curve from a P³ point:
magma Genus2Curve.m          // defines Genus2Curve; call Genus2Curve([1,2,3,4]);

// Smoke test:
magma test.m
```

Note `test.m` calls `Genus2Curve(x)` — load `Genus2Curve.m` first (or attach it)
so the function is in scope.

## Conventions & gotchas

- `τ` is stored throughout as the flat triple `[a, b, c]` representing the
  symmetric matrix `[[a, b], [b, c]]`; `SymmetricMatrix` rebuilds the matrix at
  the end.
- The surface is cut out **set-theoretically** by 3 quintics, but the saturated
  ideal's minimal basis has **more** than 3 generators (extra higher-degree
  elements). Code that needs the quintics filters with `Degree(b) eq 5` and
  asserts exactly 3 — do not assume `#MinimalBasis == 3`.
- Many sampling constants (random boxes for `zs` and `tau0`, restart counts,
  perturbation `eps`) are hard-coded — see the `TODO = make all the constants
  into parameters` markers. Convergence assumes the true `τ` falls in the
  hard-coded box.
- Precision: numeric inversion is tuned around `prec = 100..200`.
  `ReconstructCurve` needs high precision to recognize rationals; `Prec := 200`
  is the working default.
- For speed, the surface construction can be done over `GF(p)` instead of
  `Rationals()` (noted in the source) — the period computation itself is the
  expensive transcendental part.
