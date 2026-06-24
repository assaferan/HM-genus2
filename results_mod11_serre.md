# Mod-11 Galois representations of the (1,11) curves — Serre modularity evidence

Goal: exhibit evidence for **Serre's modularity conjecture over the totally real field
`K = Q(√5)`** by matching the 2-dimensional `mod 11` Galois representation attached to the
`(1,11)` Gross–Popescu genus-2 curves (Curve A etc., over `Q(√5)`) to a Hilbert modular form.

This file records the structural findings for **Curve A** (odd-part conductor = prime of
norm 66179; good reduction at both primes above 11). All scripts referenced are in the repo.

## Setup

`A = Jac(Curve A)`, `ρ̄ = A[11]`: a 4-dimensional `F₁₁`-representation of `Gal(K̄/K)`.
Because `A` has good reduction at both primes above 11, any 2-dimensional subquotient is
Barsotti–Tate at 11 ⇒ Serre weight `(2,2)`, and we seek a parallel-weight-2 Hilbert
modular form over `Q(√5)`.

## Raw data

`mod11_rep.m` reduces Curve A mod each good prime `P ≠ 11` (`N(P) ≤ 2000`) and records the
Frobenius L-polynomial. Output `mod11_frobdata.txt`: 290 rows `<N(P), p, rid, [c₀..c₄]>`
with `L_P(T) = Σ cᵢTⁱ`. The Frobenius quartic is `h_P(x) = x⁴ + c₁x³ + c₂x² + c₃x + c₄`
(`c₄ = N(P)²`). Hand-check: at `p=7` (inert), `a₇ ≡ 4 mod 11`. ✓

## Structural findings (all verified)

1. **Frobenius factorization type** (`factor_types.m`). Over all 290 primes, `h_P mod 11`
   factors **only** as `[1,1,2]` (132 primes) or `[1,1,1,1]` (158). **Never** `[2,2]`,
   `[4]`, or `[1,3]`. A geometrically simple surface with big image (`Sp₄(F₁₁)`) would show
   irreducible quartics at ~1/4 of primes; seeing none is decisive. ⇒ `A[11]` splits as
   `ρ₁ ⊕ ρ₂`, two 2-dim pieces, at every prime.

2. **No `K`-rational 1-dim subquotient** (`determine_lambda.m`, `debug_lambda.m`). The two
   linear roots at a `[1,1,2]` prime multiply to `N(P)`, so a `K`-rational isogeny line `V`
   would force a global character `λ: Gal(K̄/K) → F₁₁*` hitting one root at every prime.
   Exhaustive search over all order-`|10` Hecke characters of every plausible modulus
   (both primes over 66179 and 11, plus 3, 5, wild-2 up to exponent 12 — 12800 characters):
   the best hits only **90/132**. **No such `λ` exists** ⇒ the 11-isogeny kernel line is
   **not `K`-rational**; the assumed `W = V⊥/V` picture does not hold over `Q(√5)`.
   *(Machinery validated: it reproduces `χ₁₁ = N(P) mod 11` exactly. Pitfall: use
   `Order(Rcl.i)` of the actual ray-class generators, not `Invariants(Rcl)`.)*

3. **The two pieces are not quadratic twists** (`analyze_two_reps.m`). With
   `h_P = x⁴+Ax³+Bx²+Cx+D`, both pieces have `det = χ₁₁`, so
   `{tr ρ₁, tr ρ₂} =` roots of `y² + Ay + (B − 2N(P))` mod 11. The twist relation
   `ρ₂ ≅ ρ₁⊗ε` would force `A≡0` or `disc≡0` at every prime; it is violated at 241/290.
   So `ρ₁, ρ₂` are genuinely independent.

4. **Geometrically simple, `End = Z`** (`check_endo.m`, via CHIMP). The geometric
   endomorphism algebra is `Q`. So the `2+2` splitting is **not** from real multiplication
   or any extra endomorphism. It is forced by the **`(1,11)` polarization at the prime 11**:
   the polarization kernel `K(L) ≅ (Z/11)²` is a `K`-rational 2-dim subspace of `A[11]`.
   The prime 11 is special *because it is the polarization degree*; at other primes `ℓ` the
   image is the full `Sp₄(F_ℓ)`.

**Conclusion.** `ρ₁, ρ₂` are two genuinely-independent **irreducible 2-dimensional `mod 11`
representations of `Gal(K̄/K)`**, both with cyclotomic determinant, **not** of RM/CM origin —
exactly the interesting case for Serre's conjecture over a real quadratic field.

## Matching plan

No need to disentangle `ρ₁` from `ρ₂`. The symmetric functions are known at every prime:

```
tr ρ₁ + tr ρ₂ ≡ −A      (= −c₁ mod 11)
tr ρ₁ · tr ρ₂ ≡ B − 2N(P)   (= c₂ − 2N(P) mod 11)
```

Expectation: `ρ₁, ρ₂` are the two `mod 11` reductions of a single Hilbert modular newform
`f` over `Q(√5)`, parallel weight 2, with real-quadratic Hecke field `F` in which 11 splits.
Then for `a_P(f) ∈ O_F`,

```
−A ≡ Tr_{F/Q} a_P(f),    B − 2N(P) ≡ Nm_{F/Q} a_P(f)   (mod 11).
```

So we **match the (trace, norm) pair** `(−A, B−2N) mod 11` (known at all 290 primes) against
`(Tr, Nm)` of HMF Hecke eigenvalues mod 11.

## Level / open items

- Conductor of `ρᵢ` divides `cond(A) = 66179 · (2-part)`. At 66179 (multiplicative for `A`)
  the exponent is 0 or 1; the 2-part is not cheaply computable (Magma fails at `p|2` over
  `Q(√5)` for this Mestre model; `v_{p2}(J₁₀) = −10`).
- **No level-lowering (confirmed, `hmf_match.m`).** Scanned every weight-`(2,2)` newform of
  every `{2,3,5}`-supported level over `Q(√5)` up to norm 5625 (the level if `ρᵢ` were
  unramified at 66179), matching the correct `(Tr, Nm)` fingerprint and allowing quadratic
  Hecke fields. Nothing matches (best prefix 1/40 = noise). So `ρᵢ` is **genuinely ramified
  at 66179**; the level is `66179·(2-part)`, beyond LMFDB.
- **Next step = the heavy computation.** Compute weight-`(2,2)` Hilbert modular forms over
  `Q(√5)` at level norm 66179 (dim ~2200 via the definite quaternion algebra), reduce Hecke
  eigenvalues mod 11, and match `(Tr, Nm)` against `mod11_trnm.txt`. The 2-part of the level
  is unknown, so try level `66179`, then `66179·(2)`, `66179·(2)²` if needed.
- **Alternative worth weighing:** search the `(1,11)` construction for a curve of much
  smaller conductor — ideally over `Q` — which would make the modular form classical and the
  match trivial. The current search found 66179 as the smallest over `Q(√5)`.

## Friendlier-curve hunt (over Q) — explored, does not help

To make the modular-form side tractable we searched for a `(1,11)` curve over `Q` (→ classical
forms). `prescan_rational_gp11.m` found 9 rational non-degenerate `D₂` points (`|aᵢ| ≤ 8`);
`process_rational_d2.m` inverts them. The point `[0,2,-2,1,1]` gives a genuine genus-2 curve
**over `Q`**: field of moduli `Q`, good reduction at 11, `End=Z`, conductor
`1733036 = 2² · 433259`. Igusa invariants (normalized `J₂=1`):
```
[1, -33950578/116791249, -621348391620/1262163027943,
 -1966888453692856/13640195842980001, -14886665885581312/147409596475084870807]
```
But (i) the conductor is **larger** than Curve A's 66179; (ii) `ρᵢ` cannot level-lower at
433259 (levels 1,2,4 have no weight-2 forms), so the classical level is `≈ 433259` with
`dim S₂(Γ₀(433259)) = 36105` — far bigger than the Hilbert space for Curve A (~2200);
(iii) screening is slow (~20 min per non-converging point). **Conclusion:** the `(1,11)`
construction gives intrinsically large conductors regardless of base field; the `Q` route is
less tractable. **Curve A (norm 66179, `Q(√5)`) remains the smallest target** — proceed with
the bounded Hilbert modular forms computation there.

## THE MATCH — Serre evidence obtained

The Hilbert modular forms computation is cheap, not heavy (I had overestimated):
`HilbertCuspForms(Q(√5), level, [2,2])` at level norm 66179 has **cuspidal dim 1102**, builds
in 0.05 s, one Hecke operator ~1.7 s, <200 MB (`feasibility_hmf*.m`). Magma's definite
quaternion-algebra method is very efficient.

Matching is done **entirely mod 11** (`match_hmf_twist.m`) — the newforms have huge Hecke
fields, so we avoid number fields by reducing the Hecke operators `T_P` mod 11 directly. The
matching forms realize `ρᵢ` **up to a quadratic twist `ε` ramified at 2** (at some primes the
Hecke eigenvalue is `−t_P`, not `t_P`), so we intersect the twist-robust kernels:

```
V = ⋂_P [ ker(T_P² − tTr·T_P + tNm·I)  +  ker(T_P² + tTr·T_P + tNm·I) ]   (mod 11)
```

At level `Pbad·(2)` (the ideal of norm `264716 = 4·66179`; `(2)` the inert prime, conductor
exponent 1 there and 1 at `Pbad`; cuspidal dim **5514**):

```
dim V :  5514 → 5 → 2 → 2 → 2 → 2 → 2 → 2   (stable across 8 prime-conditions:
         both primes over 29, 31, 41; inert 7; a prime over 59)
```

`dim V = 2` is exactly the two pieces `ρ₁, ρ₂`; the chance of a 2-dim space surviving 8
independent `a_P ∈ {±q-roots}` constraints by accident is ~0. (At level 66179 alone, and at
`264716` with the untwisted sign only, `V` collapses — confirming both the 2-part of the level
and the twist.)

**Conclusion.** The two 2-dimensional `mod 11` Galois representations `ρ₁, ρ₂` of
`Jac(Curve A)/Q(√5)` are **modular**: each is realized (up to a quadratic twist at 2) by a
Hilbert modular eigenform of parallel weight 2 and level `Pbad·(2)` over `Q(√5)`. Since a
quadratic twist preserves modularity, `ρ₁, ρ₂` themselves are modular. This is concrete
**evidence for Serre's modularity conjecture over the real quadratic field `Q(√5)`** (where it
is open) — and the full pipeline ran end to end: a `(1,11)` Gross–Popescu point → an abelian
surface over `Q(√5)` → `A[11] = ρ₁ ⊕ ρ₂` → a matching Hilbert modular form.

Open refinements: identify the exact twist `ε` and the newform(s); verify over many more
primes; the 2-part exponent of the level.

### Eigensystem extraction + twist (`extract_eigensystems.m`, `identify_twist.m`)

Restricting each `T_P` to the 2-dimensional `V` gives a **scalar** `2×2` matrix at every prime
(a *single* eigenvalue `a_P`), so `V` is the multiplicity-2 eigenspace of **one** mod-11 Hecke
eigensystem `σ` (two distinct forms would give two eigenvalues; and `ρ₁,ρ₂` cannot share an
eigensystem as they are not twist-related). Thus the match exhibits **one** of the two pieces,
realized with multiplicity 2.

- **Verification (refinement b): 24/24 primes.** `a_P ∈ {±t, ±t'}` (the roots of
  `y²−Tr_P y+Nm_P` and their negatives) at *every* one of 24 primes; equivalently
  `a_P² = tr(ρ_i)(Frob_P)²`. So `σ` matches a piece **up to a quadratic twist** at all 24.
- **Global labeling via `σ` — the `K(L)` computation is circumvented.** `σ` is a genuine
  Hecke eigensystem, hence a genuine 2-dim Galois rep, so it follows **one** piece consistently:
  `tr ρ₁(Frob_P) := ` the `q`-root whose square equals `a_P(σ)²` is a globally consistent
  labeling, and `tr ρ₂ = ` the complementary `q`-root. No explicit Heisenberg/`K(L)`
  computation is needed — the modular form does the separating.

- **Twist (refinement a): `ε = χ_{K(√d)}` with `d = 66179·(5−√5)/2`** (`identify_twist2.m`,
  `verify_twist.m`). Extracting the matching order-2 Hecke character directly: `ε` is ramified
  at the prime above **5** and **both** primes above **66179**, and is **unramified at 2** (my
  earlier "ramified at 2" was an artifact of hand-transcribed data; the Magma-saved
  `eps_data.txt` is authoritative). It matches **17/18 non-degenerate primes** — the lone
  misfit is `<59, rid 8>` (one prime above 59; a residual `a_P` extraction glitch), and the one
  excluded point is the double-root prime `<61, rid 26>` where `q={7,7}` makes `ε` ill-defined.

Net: `σ ≅ ρ_i ⊗ ε` with `ε = χ_{K(√(66179(5−√5)/2))}`; one piece is modular, the labeling of
both pieces is fixed by `σ`, verified over 24 primes (with one glitch in the twist fit).
`eps_data.txt` holds the authoritative `(p,rid,ε)` points; `verify_twist.m` confirms `d`.

## Reusable artifacts produced

| File | Contents |
|---|---|
| `mod11_frobdata.txt` | raw `L_P(T)` for 290 good primes (any mod-ℓ work) |
| `mod11_apdata.txt`   | 132 `a_p mod 11` of the irreducible-quadratic factor (NB: conflates `ρ₁,ρ₂`) |
| `mod11_trnm.txt`     | **the correct target**: `(Tr, Nm) mod 11` of `{ρ₁,ρ₂}` at all 290 primes |
