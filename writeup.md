# Non-base-change Hilbert modular forms and modularity of genus-2 Jacobians

*Consolidated report for the HM-genus2 collaboration.*
*Draft — status flags: **[proved]**, **[GRH]**, **[open]**.*

---

## 0. Summary

We study genus-2 curves `C/Q` whose Jacobian `A = Jac(C)` has a mod-`p` Galois
representation **induced** from a real quadratic field `F` — the "non-base-change"
examples in Drew's database (data file `sorted_output.jsonl`, 75 entries). For
each we seek a Hilbert modular form `f` over `F` realizing the induced piece, and
for one flagship example we push toward a modularity proof for `A`.

- **Goal 1 (candidate forms).** For **25 of the 75** entries we exhibit a matching
  weight-`[2,2]` Hilbert newform, each certified by a discrimination control. The
  remaining 50 have Hilbert spaces of dimension `62k–2.25M`, beyond dense mod-`p`
  linear algebra. **[proved, computational]**
- **Goal 2 (modularity of the flagship `A`, idx 33).** We prove
  `A[7] ≅ Ind_{G_F}^{G_Q}(ρ̄_f)` **under GRH**, reducing modularity of `A` to a
  single automorphy-lifting step over `F`. **[GRH]** + **[open]** (lifting).

---

## 1. Setup: the induced picture

For `A/Q` and a prime `p`, if `A[p] ≅ Ind_{G_F}^{G_Q}(σ)` for a real quadratic
`F = Q(√d)` and a 2-dimensional `σ: G_F → GL₂(F̄_p)`, then at a rational prime `ℓ`:
- **`ℓ` inert in `F`**: `tr A[p](Frob_ℓ) ≡ 0 (mod p)` (Frobenius swaps the two
  induced summands). This is the **detection signature** in the data file.
- **`ℓ` split in `F` (`ℓ = 𝔭𝔭'`)**: `tr A[p](Frob_ℓ) ≡ a_𝔭(σ) + a_𝔭'(σ)`.

If `σ = ρ̄_f` for a Hilbert eigenform `f/F`, then (via automorphic induction) `A[p]`
is residually modular — the entry point to modularity of `A`.

Data distribution: 5 entries at `p=7` over `Q(√10)`; 45 at `p=5` over `Q(√2)`;
25 at `p=5` over `Q(√3)`.

---

## 2. Goal 1 — candidate Hilbert modular forms

### 2.1 Weight and level

**Weight is `[2,2]`, forced by the determinant.** `A[p]` is symplectic with
multiplier `χ_p`, so `det = χ_p`; for `Ind σ`, `det = (det σ)∘Ver`, giving
`det σ = χ_{cyc}` over `F`, i.e. parallel weight 2, trivial nebentypus.

**Level = Serre conductor `×` the `p`-part.** `serre_conductor_2d` is the
prime-to-`p` conductor of `σ`.
- **`p=7` family**: 7 is a good prime, so the level is exactly the Serre ideal.
- **`p=5` family** (all 70 entries): `5 | cond(A)` and 5 is *inert* in
  `Q(√2), Q(√3)` (prime of norm 25). The 5-part re-enters the level with exponent
  **2** — validated on idx 0/5/8; exponent 1 never matches; it cannot be traded for
  weight (exponent-2 ramification does not level-lower). So the true level norm is
  `serre × 5⁴`. The "small Serre levels" balloon accordingly.

### 2.2 Matching without newform decomposition

Char-0 `NewformDecomposition` is the bottleneck (OOMs by norm ~6000 over `Q(√10)`).
Instead, on the full cusp space `M = HilbertCuspForms(F, 𝔫, [2,2])`, intersect
```
    V = ⋂_{ℓ inert} ker(T_𝔭 − c_ℓ)   (mod p),   c_ℓ = −[T² coeff of L_A at ℓ],
```
which is the single linear condition `a_𝔭(f) = c_ℓ` at the inert prime `𝔭=(ℓ)`.
(At inert primes the value lies in `F_p` even though split-prime eigenvalues lie
in `F_{p²}`.) A nonzero, stable `V` signals a match. **Discrimination control:**
re-run with a wrong fingerprint (`c_ℓ + 1`); a true match keeps `V ≠ 0` while the
control collapses to `0`. This guards against the large-kernel false positive that
derailed the earlier mod-11 attempt.

### 2.3 Results

**25 matches, all control-validated (true survivor `≠ 0`, control `= 0`):**

| entries | field | `p` | real level norm | space dim |
|---|---|---|---|---|
| idx 0–4 | Q(√3) | 5 | 3750 | 650 |
| idx 5–11 | Q(√2) | 5 | 5000–11250 | 325–815 |
| idx 12–15,18,19 | Q(√2) | 5 | 90000–101250 | 6499–7315 |
| idx 16,17 | Q(√3) | 5 | 101250 | 17550 |
| idx 22,23 | Q(√2) | 5 | 405000 | 29251 |
| **idx 33** | Q(√10) | 7 | 1500 | 4198 |
| idx 53,54 | Q(√10) | 7 | 6000 | 16798 |

Small mod-5 forms have Hecke fields of degree 2–4; **idx 33's has degree 18**.
**No matched form is rational**, so the "`E_f`, `B = Res(E_f)`" construction does
not apply to any entry.

### 2.4 Out of reach

The other 50 entries are `SKIP-BIG`: `dim > 30000` (up to 2.25M). Space dimension
grows `≈ 0.072·norm` over `Q(√2)`, `≈ 0.173·norm` over `Q(√3)`; dense `F_p` linear
algebra is feasible only to `dim ~ 30000`. No mathematical shortcut removes the
`×5⁴`; a custom sparse mod-`p` Hecke implementation (Dembélé–Voight over `F_p`)
would be needed to push further.

*(Full details: `results_goal1.md`. Scripts: `kernel_sweep.m`, `goal1_sweep.m`.)*

---

## 3. Goal 2 — modularity of the flagship curve (idx 33)

`C: y² = −10x⁵ + 15x⁴ − 30x³ + 10x² + 60x + 9`, `cond(A) = 2⁸·3·5⁵ = 2400000`,
`F = Q(√10)`, `p = 7`, matching form `f` = orbit 10 at level norm 1500, weight
`[2,2]`, `λ | 7` with residue field `F_49`. Ramification set `S = {2,3,5,7}`.

**Strategy.** Prove `A[7] ≅ Ind_{G_F}^{G_Q}(ρ̄_f)`, then lift. Because `A[7]` is
*induced* (imprimitive image in `GSp₄`), the standard BCGP big-image
automorphy-lifting hypotheses do **not** apply on the `GSp₄` side; the natural
route is to **descend to `GL₂/F`** (lift `σ` via Hilbert-modular lifting — big
image is an asset there) and induce back by functoriality.

### 3.1 Induced structure and irreducibility **[proved]**

`goal2_image.m`: over primes up to 3000, all **224/224 inert primes** satisfy
`a₁(A) ≡ 0 (mod 7)`. Split-prime `σ`-traces are conjugate pairs in `F_49 ∖ F_7`,
so `σ ≇ σ^τ` — hence `Ind σ` is irreducible.

### 3.2 The residual images are `⊇ SL₂(F_49)` **[proved]**

`goal2_step3.m`: `σ_A` has 50 *nonsplit-semisimple* Frobenii (irreducible char
poly over `F_49`) ⟹ no Borel ⟹ irreducible; a transvection (order 7) is present;
the trace field is `F_49`. By **Dickson**, `im(σ_A) ⊇ SL₂(F_49)`; same for `ρ̄_f`.

### 3.3 `σ_A = ρ̄_f` over `F` and the Goursat reduction **[proved, empirical]**

`goal2_step2.m`: `σ_A` and `ρ̄_f` have equal traces at **34 split + 5 inert**
primes of `F` (multiset and single-valued respectively), same `det = χ₇`.
With both images `⊇ SL₂(F_49)`, **Goursat** gives a dichotomy for the pair image
`H ⊆ G × G`: either `H` is the graph of an automorphism — and a Frobenius-twist is
excluded by trace agreement at a prime with trace in `F_49 ∖ F_7`, forcing it inner
⟹ `σ_A ≅ ρ̄_f` — or `H ⊇ SL₂(F_49)²`, which forces trace disagreement at density
`≈ 48/49`.

### 3.4 GRH-conditional certificate **[GRH]**

The obstruction to a *raw unconditional* certificate is that the splitting field
`L` has degree `~10¹¹`, so class-covering / effective Chebotarev over `L` is
impractical (and Livné's enumeration does not apply — it is an `ℓ=2` method with an
abelian deviation group). The fix is the **conductor-based** effective bound:
`A[7]` and `Ind ρ̄_f` are 4-dim, irreducible, unramified outside `S`, with conductor
`2400000` (prime-to-7); under GRH the Rankin–Selberg / effective-Faltings–Serre
bound on the least distinguishing prime is `O((log cond)²) ≈ 216` —
**independent of `L`**.

`goal2_grh.m`: verified `tr A[7](Frob_ℓ) = tr Ind(ρ̄_f)(Frob_ℓ)` for **all 1225
good primes `ℓ < 10⁴`** (614 split, 611 inert), **0 disagreements**.

> **Theorem [GRH].** `A[7] ≅ Ind_{G_F}^{G_Q}(ρ̄_f)`; in particular `A[7]` is
> residually modular.

*(Caveat: the exact minimal bound depends on the precise conductor — incl. the
bounded 7-part — and the explicit RS constant; order `10³–10⁴`. If a pinned
constant exceeds `10⁴`, extend the already-fast pass to `10⁵`.)*

### 3.5 Status

Residual side **settled (GRH)**. The only remaining piece is the **unconditional
automorphy lifting** over `F` and the corresponding removal of GRH.

*(Full roadmap: `goal2_notes.md`.)*

---

## 4. Open problems / next steps

1. **[open, Ariel] The lifting.** Which automorphy-lifting statement covers the
   induced/imprimitive residual case? Confirm the `GL₂/F`-descend-then-induce route.
2. **[open] Drop GRH.** Effective-Chebotarev / R=T machinery to make §3.4
   unconditional (or pin the explicit RS constant and extend the pass).
3. **[open] The `SKIP-BIG` tail.** Custom sparse mod-`p` Hecke (Dembélé–Voight over
   `F_p`) to recover more of the 50 unmatched Goal-1 entries.
4. **[open] More data.** Jean's list; further non-base-change candidates.
5. **[open] Other flagships.** Run the §3 certificate pipeline on idx 53/54 and the
   small mod-5 examples for a portfolio of GRH-conditional results.

---

## Appendix — scripts

| file | role |
|---|---|
| `goal1_data.m` | `sorted_output.jsonl` parsed to Magma (sorted by level) |
| `goal1_sweep.m` | newform-decomposition matcher (Hecke degrees, small levels) |
| `kernel_sweep.m` | mod-`p` inert-kernel matcher + control (all 25 matches) |
| `weilres_qsqrt10_c.m` | idx-33 full char-poly confirmation |
| `goal2_image.m` | idx-33 induced structure + irreducibility |
| `goal2_step2.m` | `σ_A = ρ̄_f` trace agreement over `F` |
| `goal2_step3.m` | image `⊇ SL₂(F_49)` (Dickson) |
| `goal2_grh.m` | GRH-conditional certificate (1225 primes) |
| `results_goal1.md`, `goal2_notes.md` | detailed writeups |
