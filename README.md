# HM-genus2

Reconstruct a genus-2 curve from a point parametrizing an abelian surface via the
**Horrocks–Mumford (HM) bundle**, using [Magma](http://magma.maths.usyd.edu.au/).

Given `x ∈ P³` (a section of the HM bundle), the pipeline recovers a genus-2 curve
`C/ℚ` whose Jacobian is linked to the HM surface `A_x` by a **rational cyclic
5-isogeny**:

```
x ∈ P³  ──HMSurface──►   A_x ⊂ P⁴   (smooth, degree 10, (1,5)-polarized)
        ──x_to_tau──►    period matrix τ ∈ H₂   (numeric inversion of the theta map)
        ──Genus2Curve──► genus-2 curve C/ℚ  with a rational 5-isogeny on Jac(C)
```

## Mathematical sketch

- The **Horrocks–Mumford bundle** `F` on `P⁴` has `h⁰(F)=4`; a section
  `s = Σ xᵢ·secᵢ` (`x ∈ P³`) cuts out a `(1,5)`-polarized abelian surface
  `A_x ⊂ P⁴` of degree 10.
- `A_τ = ℂ²/(τℤ² + diag(1,5)ℤ²)` embeds in `P⁴` by the five `(1,5)` theta
  functions. The forward map is transcendental with no closed form, so `x → τ`
  is solved **numerically**: a point `z ∈ A_τ` maps to a point that must satisfy
  the defining quintics of `A_x`, giving a residual `q(Θ(z))` minimized over the
  three entries of `τ`.
- The principally polarized `[τ | I]` is `5`-isogenous to `A_x`; its Jacobian is
  the reconstructed curve. The curve is recovered from its (rational) **absolute
  Igusa invariants**.

## Files

| File              | Role |
|-------------------|------|
| `HMsurface.m`     | HM bundle and the surface `A_x = HMSurface(F, x)`; the `(1,5)` theta embedding. |
| `Inversion.m`     | Reference inverse `x_to_tau` (finite-difference Gauss–Newton) + `verify_tau`. |
| `InversionFast.m` | Production inverse `x_to_tau_fast`: analytic Jacobian, Levenberg–Marquardt, adaptive truncation, two-phase (low→high) precision with stagnation/eigenvalue pruning. |
| `Genus2Curve.m`   | Top-level driver `Genus2Curve(x : K := ...)`: `x → τ → curve/K` via Igusa invariants recognized in `K` (default `K = ℚ`; any number field supported). |
| `Reduction.m`     | Good-reduction-at-5 filter, conductor (per-prime over a number field), minimal quadratic twist, Frobenius polynomials. |
| `Heights.m`       | Enumerate `P³(K)` points (`K` quadratic) by increasing Weil height. |
| `search_quadratic.m` | Search driver: screen `ℚ(√2)`-points by increasing height for good reduction at 5, rank by conductor size. |
| `test.m`          | Checks the reconstructed curve has a **rational cyclic 5-isogeny**. |
| `test_quadratic.m`| Example: reconstructs a curve over `ℚ(√2)` from `x = [1, √2, 3, 4]`. |

## Dependencies

- **Magma** (commercial CAS) — provides the HM bundle, theta sums, scheme/ideal
  machinery, and curve/Jacobian arithmetic.
- **[CHIMP](https://github.com/edgarcosta/CHIMP)** — must be checked out as a
  sibling directory `../CHIMP` (so that `AttachSpec("../CHIMP/CHIMP.spec")`
  resolves). Provides `ComplexFieldExtra`, `ReduceSmallPeriodMatrix`,
  `ThetaDerivatives`, etc. CHIMP is **not** redistributed here.

Expected layout:

```
parent/
├── HM-genus2/      ← this repo
└── CHIMP/          ← https://github.com/edgarcosta/CHIMP
```

## Usage

```magma
// Build the curve from a P^3 point:
magma
> load "Genus2Curve.m";
> C := Genus2Curve([1,2,3,4]);
> C;
Hyperelliptic Curve defined by y^2 = 12*x^6 + 4*x^5 + 36*x^4 + 12*x^3 + 39*x^2 + 8*x + 16 over Rational Field

// Run the 5-isogeny test:
magma Genus2Curve.m test.m
```

### Over a number field

The same driver handles `x ∈ P³(K)` for a number field `K`: pass `K := <field>` and
`x` over `K`. The surface is built over `K`, inverted via a chosen complex embedding,
the Igusa invariants are recognized in `K`, and the curve is descended to `K`
(Mestre's algorithm):

```magma
> load "Genus2Curve.m";
> K<s2> := QuadraticField(2);
> C, QI := Genus2Curve([K| 1, s2, 3, 4] : K := K);   // Embedding := 1 -> s2 |-> +1.41421...
> BaseRing(C);     // Q(sqrt2) -- the curve descends (no Mestre obstruction here)
```

or `magma Genus2Curve.m test_quadratic.m`. Notes:
- Over a number field only some embeddings yield a positive-definite period τ (which
  ones is `x`-dependent). `Embedding := 0` (default) auto-tries the embeddings and
  uses the first that works; `Embedding := i` forces the `i`-th (roots of `K`'s
  defining polynomial, sorted canonically). Different working embeddings give
  Galois-conjugate curves.
- The descended model can have very large coefficients — Magma's genus-2 model
  reduction does not support number fields, so reduction is applied only when `K = ℚ`.
- The Igusa invariants always lie in `K` (the field of moduli); a *model* over `K`
  exists exactly when the Mestre conic has a `K`-rational point.

### Good reduction at 5 and conductor

`Reduction.m` provides:

```magma
import "Reduction.m": PotentialGoodReductionAt5, GoodReductionAt5, ReductionReport;
PotentialGoodReductionAt5(QI, K);   // cheap, model-free filter from the Igusa invariants
GoodReductionAt5(C);                // definitive (conductor exponent 0 at p | 5)
ReductionReport(C);                 // status at 5 + conductor exponents per bad prime
```

- `PotentialGoodReductionAt5` uses the Igusa valuation criterion
  (`v_p(J₁₀)/10 ≤ v_p(J_{2i})/2i`) on the invariants the pipeline already returns —
  so you can **screen many `x` for good reduction at 5 without building the curve model**.
- The actual conductor exponent comes from `Conductor(C, p)`, which over a number field
  works at all **odd** primes (in particular `p | 5`) but is **not implemented at
  `p | 2`** ("fibre blowups"). Over `ℚ` the full `Conductor(C)` is available (its
  2-part uses Ogg's formula heuristically when `v₂(disc) ≥ 12`).
- Note: by construction every curve from this family already has a rational 5-isogeny,
  so the search is really over **good reduction at 5 + small conductor**.

### Searching a family

`search_quadratic.m` enumerates `ℚ(√2)`-points of `P³` in **increasing height order**
(`Heights.m`: `BoundedHeightPoints`), and for each computes the Igusa invariants
(`Genus2Invariants` — the period computation **without** the Mestre descent, the cheap
way to screen), keeps only those with good reduction at 5, and ranks survivors by a
conductor-size proxy (norm of the bad primes away from 2 and 5, found cheaply at small
primes). Results stream to `search_results.txt`. For the best survivors, build the model
with `Genus2Curve` and read exact conductor exponents (odd primes) via `ReductionReport`.

The enumeration **excludes points `x ∈ P³(ℚ)`** (`ExcludeRational`), since a rational `x`
makes the whole construction — and the curve — defined over ℚ (so `C/K` is a base change).
As a complete check the driver also **rejects curves with rational Igusa invariants**
(geometrically defined over ℚ), which catches the rare irrational `x` landing on a ℚ-curve.

For each point the driver: computes `τ`; filters by **potential good reduction at 5**
(invariant criterion) and by **`End(J)=ℤ`** (`GeometricEndomorphismDimension(τ)==1`,
i.e. geometrically simple — both cheap and twist-invariant); and for survivors descends
to the curve, applies the **minimal quadratic twist** (`MinimalTwist`: twist to good
reduction at every prime of potential good reduction, so the conductor is as small as
possible), and records the conductor exponents and the factored **Frobenius polynomial
at 5**.

Low height tends to give low conductor, so processing in height order surfaces the best
candidates first. **Scale caveat:** `P³(K)` has ~`B⁸` points of height ≤ `B`, so an
exhaustive height bound is infeasible beyond tiny `B`; the driver enumerates within a
coordinate box (`BoxM`) and processes the lowest `MaxPts` points by height (both
configurable). The per-point period computation is the cost; the filters are cheap and
run before the (expensive) Mestre descent. Many of the very lowest-height points are
degenerate (the section leaves the abelian locus) and are skipped automatically.

**Twist subtleties.** The reconstructed curve is only determined up to quadratic twist,
and good reduction at 5 / the conductor / the Frobenius polynomial at 5 are all
twist-dependent (the 5-isogeny structure and `End=ℤ` are twist-invariant). `MinimalTwist`
pins the model with smallest conductor. Worked example: `x = [0,-1,1,1]` (height 1) gives
a **simple** Jacobian with a rational 5-isogeny whose minimal model has **good reduction
away from 2** (in particular good at 5), with Frobenius polynomial at 5 equal to
`(25T²+4T+1)²` (so its reduction mod 5 is `≅ E²` over `𝔽₂₅`).

## Notes

- The numeric inversion uses random restarts (only ~10% of starts land in the
  true basin), so a run does several short low-precision searches before a single
  full-precision polish. `Genus2Curve(x : Prec := 300)` controls the working
  precision.
- The test checks a **rational 5-isogeny** (a Galois-stable order-5 subgroup of
  `J[5]`), *not* a rational 5-torsion point — the latter is strictly stronger and
  need not hold. Concretely: at every prime `p` of good reduction (`p ≠ 5`), the
  mod-5 reduction of the `L`-polynomial of `C/𝔽ₚ` must have a root in `𝔽₅`.

## License

MIT — see [LICENSE](LICENSE).
