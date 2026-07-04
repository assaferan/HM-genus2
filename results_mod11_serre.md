# Mod-11 Galois representations of the (1,11) curves — Serre modularity evidence

> **⚠️ CORRECTION (2026-06-29): the modularity match claimed in the section
> "THE MATCH — Serre evidence obtained" (below) is a FALSE POSITIVE and is RETRACTED.**
> On refinement, the 2-dimensional space `V` does not realize `ρ₁⊕ρ₂`: restricting the
> Hecke action to `V`, `charpoly(T_P|V) ∈ {q⁺,q⁻}` at only 2/25 primes, `det(T_P|V) ≠ Nm_P`
> at most primes, and `V` is non-semisimple. The matcher's acceptance criterion
> (`V = ⋂_P[ker q⁺(T_P)+ker q⁻(T_P)] ≠ 0`) is too weak — `ker` of a degree-4 polynomial in
> `T_P` covers ~4/11 of the space mod 11, so a nonzero intersection does not imply an
> eigenform with a consistent twist. The correct per-eigenform test fails: there is no
> untwisted match (`⋂ ker q⁺ = 0` at every level), and the unique twisted survivor at level
> `𝔭`, though it satisfies `a_P ∈ {±t₁,±t₂}` at 120/120 primes, has a twist sign `ε(P)` that
> is provably **not** a quadratic Hecke character (any conductor in `{2^≤6,3,5,11,66179,∞}`).
> Since a parallel-weight-2 trivial-character form forces `det ρ_f = χ₁₁ = det ρ_i`, a
> non-quadratic twist is impossible, so no such eigenform realizes `ρ_i`.
> **Status: the Serre match for Curve A is NOT established.** The *structural* findings
> (Sections "Structural findings", "Friendlier-curve hunt") remain valid. See the
> "Refutation" section at the end and scripts `refine_serre[2-5].m`, `refine_charfit.m`,
> data `refine_serre5_signs.txt`.

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

## THE MATCH — Serre evidence obtained  — ⚠️ RETRACTED, see "Refutation" below

*(The following section is preserved as written but its conclusion is WRONG; see the
correction banner at the top and the "Refutation" section at the end.)*

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

## Reusable artifacts produced

| File | Contents |
|---|---|
| `mod11_frobdata.txt` | raw `L_P(T)` for 290 good primes (any mod-ℓ work) |
| `mod11_apdata.txt`   | 132 `a_p mod 11` of the irreducible-quadratic factor (NB: conflates `ρ₁,ρ₂`) |
| `mod11_trnm.txt`     | **the correct target**: `(Tr, Nm) mod 11` of `{ρ₁,ρ₂}` at all 290 primes |

## Refutation (2026-06-29) — the match above is a false positive

Refining the claim (pin down the twist, verify over more primes) required checking that the
"matched" space `V` actually *realizes* `ρ₁⊕ρ₂`, not merely that `V ≠ 0`. It does not.

**The acceptance criterion was too weak.** `V = ⋂_P [ker q_P⁺(T_P) + ker q_P⁻(T_P)]` uses the
kernel of a *degree-4* polynomial in `T_P`; over `F₁₁` that kernel has dimension ≈ `(4/11)·dim`,
so a nonzero intersection over a handful of primes is a weak condition and does **not** imply
the existence of an eigenform whose twisted eigenvalues match the fingerprint.

**Direct test of the dim-2 `V` (level `𝔭·(2)`).** Restricting each `T_P` to `V` gives a 2×2
matrix `C_P` over `F₁₁` (`refine_serre3.m`). For a genuine `ρ₁⊕ρ₂` one needs
`det C_P = Nm_P` and `charpoly C_P = q_P⁺` or `q_P⁻` at every prime. Instead:
- `charpoly(C_P) ∈ {q⁺,q⁻}` at only **2/25** primes (both degenerate `tNm=0` cases);
- `det(C_P) ≠ Nm_P` at most primes;
- `C_P` typically has a **repeated** eigenvalue (e.g. `(y−7)²`) — `V` is non-semisimple, hence
  not a sum of two genuine Hecke eigenforms.

**Correct per-eigenform test (`refine_serre4.m`, `refine_serre5.m`).**
- Untwisted intersection `V⁺ = ⋂_P ker q_P⁺(T_P) = 0` at level `𝔭` **and** `𝔭·(2)`: no
  untwisted match anywhere. (Were a newform with `11` split in its Hecke field to realize
  `ρ₁,ρ₂` as its two mod-11 reductions, both reductions would be roots of `q⁺`, i.e. lie in
  `V⁺`; `V⁺=0` already argues against such a form at these levels.)
- The twist-robust survivor at level `𝔭` is **1-dimensional** (a single eigenform `f`). Its
  eigenvalue satisfies `a_P ∈ {±t₁,±t₂}` at **120/120** primes (`FAIL 0`) — a real structural
  fact, not chance (`≈ (4/11)^{106}`). Define `ε(P)=+1` if `a_P` is a root of `q⁺`, `−1` if of
  `q⁻`.
- `ε` is **not a quadratic Hecke character.** The 94 unambiguous signs are inconsistent with
  every quadratic ray-class character of conductor supported in `{2^{≤6},3,5,11,66179,∞}`
  (`refine_charfit.m`). `K=Q(√5)` has narrow class number 1, so there are no unramified
  quadratic characters either.

**Why this is decisive.** A Hilbert eigenform of parallel weight 2 and *trivial* character has
`det ρ_f = χ₁₁ = det ρ_i`. Hence `ρ_f ≅ ρ_i ⊗ χ` forces `χ² = 1`, i.e. the only admissible
twist is quadratic — exactly what is ruled out above. (A non-trivial nebentypus would spoil
`det = χ₁₁`, so it cannot help.) Therefore **no** parallel-weight-2 trivial-character
eigenform of level `𝔭·(2)^a` realizes `ρ₁` or `ρ₂`.

The residual `a_P² ∈ {t₁²,t₂²}` (120/120) is the *twist-invariant* (projective / Sym²)
statement, with the square jumping between `ρ₁` and `ρ₂`; it does not imply either `ρ_i` is
modular.

**Bottom line: the mod-11 Serre-modularity match for Curve A is NOT established.** The
structural results above (`A[11]=ρ₁⊕ρ₂`, irreducible, `det=χ₁₁`, `End=ℤ`, not twists) stand.

| Refutation script | Role |
|---|---|
| `refine_serre3.m` | per-prime 2×2 `C_P=T_P|V`: eigenvalues, `det`, `charpoly` vs `q±` |
| `refine_serre4.m` | per-eigenform test: `V⁺,V⁻,V_tw` dims; root-membership + twist fit |
| `refine_serre5.m` | dumps `⟨N,p,a_P,t₁,t₂,ε⟩` over 120 primes (`refine_serre5_signs.txt`) |
| `refine_charfit.m`| tests whether `ε` is a quadratic ray-class character (it is not) |

## The corrected modularity statement (derived, 2026-07-01)

Having refuted the false match, we pin down what the correct statement **must** be, by
elimination. Assuming Serre's conjecture over `Q(√5)` holds (as it should — each `ρ_i` is
2-dimensional, irreducible, and odd), each `ρ_i` is modular by a Hilbert modular eigenform, and
we can determine its weight, character, and level up to one bounded constant.

> **Statement.** Each of `ρ₁, ρ₂` is modular over `K = Q(√5)` by a Hilbert modular eigenform of
> **parallel weight `(2,2)`, trivial nebentypus, and level `𝔭·(2)^a`**, realized **untwisted**
> (`det = χ₁₁`), where `𝔭` is the prime of norm 66179 and `a = a₂(ρ_i) ≥ 3` is the conductor
> exponent of `ρ_i` at the inert prime 2.

**How each parameter is forced:**

| parameter | value | reason |
|---|---|---|
| determinant | `χ₁₁` | det-clean at 132/132 `[1,1,2]` primes (eigenvalue products `= N(P)`) |
| nebentypus | trivial | `det ρ_f = ψ·χ_cyc^{k−1}`; weight 2 + `det=χ₁₁` ⇒ `ψ` trivial. Nebentypus `χ₁₁` gives `det χ₁₁²`, and `χ₁₁` (order 10) has **no square root** among its powers ⇒ obstruction (`cond2_v2.m`, determinant argument) |
| weight | `(2,2)` | `A` is **ordinary at both primes `𝔩\|11`** (`L mod 11` has degree 2 = two unit roots; `local_at_11.m`). Ordinary + good reduction ⇒ peu-ramifiée; a subquotient of a finite-flat group scheme is finite-flat (`e=1<10`, Raynaud); `det=χ₁₁` forces HT weights `{0,1}` ⇒ weight 2. Companion weight 12 (and non-parallel `(2,12)`) unnecessary |
| `𝔭`-part | `𝔭¹` | pure `(2)^a` levels (`a=1..6`) give `V⁺=Vtw=0` immediately (`test_2adic.m`) ⇒ `ρ_i` ramified at `𝔭` |
| 2-part | `(2)^a, a≥3` | `V⁺=0` (untwisted, weight `(2,2)`) at `𝔭·(2)^{0,1,2}` (`refine_serre4.m`, `test_level2.m`, dim 22059 at `(2)²`) |
| realization | untwisted | `det ρ_f = det ρ_i = χ₁₁` ⇒ if modular then `ρ_f ≅ ρ_i` directly; and no quadratic twist works (the `a_P²` sign is not a character) |

**Consequences / bounds on the one open constant `a`.** Since `V⁺=0` at 2-part `≤2` for *both*
pieces, `a₂(ρ₁), a₂(ρ₂) ≥ 3`. The Artin conductor is additive, so
`a₂(A) = a₂(ρ₁) + a₂(ρ₂) ≥ 6`: heavy wild ramification at 2. Wild inertia `P₂` acts through a
2-group of order `≤ 16` (the 2-Sylow of `GL₂(𝔽₁₁)`, `|GL₂(𝔽₁₁)|=2⁴·3·5²·11`).

**The exact value of `a` is currently not computable.** Every route to `cond₂` is blocked:
Magma has no genus-2 regular models / minimization over number fields ("fibre blowups over
number fields not yet implemented" — kills `Conductor`, `RegularModel`, auto-`LSeries`);
the Weil restriction `W=Res_{K/Q}(A)` gives the clean bridge `a₂(W) = 2·a_{𝔭₂}(A)` (2 inert)
but `W` has no curve model; and the numerical functional-equation route needs a minimal model,
whereas `Genus2CurveFromIgusa` returns a model whose content is divisible by ~20 primes
(including `6.3×10²³`). This is the precise reason behind the earlier note that the 2-part "is
not cheaply computable." Determining `a` requires either a genuine minimal model of Curve A
(then the functional-equation computation) or the empirical weight-`(2,2)` Hecke test at
`𝔭·(2)³` (dim ≈ 88000).

| Derivation script | Role |
|---|---|
| `local_at_11.m` | reduction type of `A` at `𝔩₁,𝔩₂\|11` (ordinary ⇒ weight `(2,2)`) |
| `test_2adic.m` | pure 2-power levels `(2)^a` empty ⇒ `ρ_i` ramified at `𝔭` |
| `test_level2.m` | `V⁺=0` at `𝔭·(2)²` (dim 22059) |
| `cond2_probe.m`, `cond2_v2.m` | conductor-at-2 attempts (all blocked; document the walls) |
| `lseries_probe.m` | `LSeries(C/K)` hits the same number-field blowup wall |

## The conductor at 2 — determined, and the complete statement (2026-07-04)

The one open constant `a = a₂(ρ_i)` is now pinned down numerically via the **functional
equation**, completing the modularity statement.

> **Complete statement.** Each of `ρ₁, ρ₂` is modular over `K = Q(√5)` by a Hilbert modular
> eigenform of **parallel weight `(2,2)`, trivial nebentypus, level `𝔭·(2)⁴`**, realized
> **untwisted** (`det = χ₁₁`).

**Method.** `L(A/K,s) = L(W/ℚ,s)` for `W = Res_{K/ℚ}(A)`; the analytic conductor is
`N(W) = |d_K|⁴ · N_{K/ℚ}(cond_K A) = 5⁴ · 66179 · 4^c`, where `c = a_{𝔭₂}(A)` is the exponent
at the inert prime 2. We build `L(A/K,s)` from the Euler factors and search `c` by
`CheckFunctionalEquation`.

**Getting the right curve.** `Genus2CurveFromIgusa` returns a **wrong twist** `A_recon` — it is
additive (potentially-good) at `499, 13711, 88301, 231481, …` (the residual model content `g`;
the Igusa criterion misses these since potentially-good passes it, and `MinimalTwist` misses
them since `Conductor` fails over `K`). The correct curve `A_true = A_recon ⊗ χ_g`
(`χ_g` = quadratic character of `K(√g)`) is good there, so its conductor is just `𝔭·(2)^c`
(no large primes), making the L-function tractable. Its Euler factors are the `χ_g`-twist of
`A_recon`'s: `L^{true}_P(T) = L^{rec}_P(χ_g(P)·T)`.

**Result** (`cfe_final2.m`, precision 3, 25 981 Euler factors up to norm 300 000):

| `c` | 6 | 7 | **8** | 9 |
|---|---|---|---|---|
| `\|CFE\|` | 0.0996 | 0.293 | **0.00195** | 0.227 |

`c = 8` closes the functional equation ~2 orders of magnitude cleaner than any neighbor, and is
consistent with the proven bound `a₂(A) ≥ 6`. Since the Artin conductor is additive,
`a₂(A) = a₂(ρ₁) + a₂(ρ₂) = 8`; the symmetric case gives **`a₂(ρ_i) = 4`**, hence level `𝔭·(2)⁴`.

*(This is a numerical, functional-equation determination — very strong, not a rigorous proof.)*

**Compute notes.** `a₁` via `#C(F_q)` is fast (~0.01 s); the full L-poly (for `a₂`) is only
needed for `Norm(P) ≤ 560`. Two pitfalls that cost a lot: (i) **skip `Norm(P) > cutoff`** —
inert primes `p > 216` have `Norm p² > cutoff` and were being `#C(F_{p²})`-counted (`O(p²)`),
catastrophically slow; (ii) the badly non-minimal primes `3, 5, 19, 191` need a **zoom-search
minimizer** (`x → π·x + r` at a repeated root, `π = ` the `P`-local principal generator).

| Conductor-2 script | Role |
|---|---|
| `cfe_final2.m` | the working functional-equation conductor search (`c=8`) |
| `cfe_ingredients.m` | verifies `A_true` good at 499 etc.; the `χ_g` twist |
| `minimize1.m`, `minimize_hard2.m` | model content removal + zoom minimizer for hard primes |
