# Goal 1 — candidate Hilbert modular forms for the non-base-change list

**Input:** `sorted_output.jsonl` (Ariel's data file), 75 genus-2 curves from
Drew's database flagged as likely giving **non-base-change** HMFs — i.e. whose
mod-`p` Galois representation `A[p]` is *induced* from a real quadratic field
`F = Q(√d)`. Detection signature (in the file's `traces`): `a_ℓ(A) ≡ 0 (mod p)`
at every prime `ℓ` inert in `F` (Frobenius at an inert prime swaps the two
induced summands, so the trace vanishes).

**Task:** for each entry find a Hilbert eigenform `f` over `F`, weight `[2,2]`,
whose mod-`p` reduction realizes the 2-dimensional piece, i.e.
`A[p] ⊗ F̄_p ≅ Ind_{G_F}^{G_Q} ρ̄_{f,λ}`.

Distribution: 5 entries mod 7 over `Q(√10)`; 45 mod 5 over `Q(√2)`; 25 mod 5
over `Q(√3)`.

---

## Results: 25 / 75 matched (all validated), 50 out of computational reach

Every match below was confirmed by a **discrimination control**: on the same
Hecke matrices, the true fingerprint `{c_ℓ}` yields a **nonzero** surviving
mod-`p` eigenspace while a deliberately wrong fingerprint (`c_ℓ+1`) collapses to
**0**. This rules out the oldform/large-kernel false positive that derailed the
mod-11 Curve-A attempt.

| idx | field | p | real level norm | space dim | survivor (true/ctrl) | Hecke deg† |
|----:|-------|--:|----------------:|----------:|:--------------------:|:----------:|
| 0–4  | Q(√3) | 5 | 3750    | 650   | 20 / 0 | 4 |
| 5–7  | Q(√2) | 5 | 5000    | 325   | 12 / 0 | 2 |
| 8–10 | Q(√2) | 5 | 10000   | 649   | 36 / 0 | 2 |
| 11   | Q(√2) | 5 | 11250   | 815   | 12 / 0 | 4 |
| 12   | Q(√2) | 5 | 90000   | 6499  | 60 / 0 | — |
| 13–15,18,19 | Q(√2) | 5 | 101250 | 7315 | 24–40 / 0 | — |
| 16,17 | Q(√3) | 5 | 101250 | 17550 | 24 / 0 | — |
| 22,23 | Q(√2) | 5 | 405000 | 29251 | 48 / 0 | — |
| **33** | Q(√10) | 7 | 1500  | 4198 | 8 / 0 | **18** |
| 53,54 | Q(√10) | 7 | 6000  | 16798 | 32 / 0 | — |

† Hecke-field degree, where the newform decomposition reached (small levels).
Blank = found by the kernel method, which does not compute the Hecke field.

**idx 33 is the flagship** (Ariel's `y² = −10x⁵+15x⁴−30x³+10x²+60x+9`): beyond
the trace match, the **full degree-4 characteristic polynomial mod 7** of `A`
agrees with `Ind ρ̄_f` at every good prime tested up to ~110 (split *and* inert,
0 mismatches). The matching form is **orbit 10, Hecke field degree 18**, reduced
at `λ|7` with residue field `F_49` — **not rational**.

The 50 unmatched entries are `SKIP-BIG`: the relevant Hilbert space has dimension
`> 30000` (up to 2.25 million), where the dense `F_p` linear algebra is
infeasible. See feasibility note below.

---

## Method

**1. Weight is `[2,2]`, forced by the determinant.** `A[p]` is symplectic with
multiplier `χ_p`, so `det = χ_p`. For `Ind σ`, `det = (det σ)∘Ver`, so
`det σ = χ_{cyc}` over `F` ⟹ parallel weight 2, trivial nebentypus.

**2. Level = Serre conductor × the p-part.** The tabulated `serre_conductor_2d`
is the conductor of `σ` **prime to p** (Serre convention).
- For the **mod-7** family, `7 ∤ cond(A)` (7 is a good prime), so the level is
  the Serre ideal exactly.
- For **every mod-5** entry, `5 | cond(A)`, and 5 is **inert** in `Q(√2)`,`Q(√3)`
  (prime of norm 25). The 5-part must be restored into the level with exponent
  **2** (validated on idx 0/5/8; exponent 1 never matches). It cannot be traded
  for weight, because exponent-2 ramification does not level-lower. Hence the
  true level norm = `serre × 625` — the "small Serre levels" balloon.

**3. Matching without newform decomposition (the enabling step).** Char-0
`NewformDecomposition` is the bottleneck (it OOM'd already at norm ~6000 over
`Q(√10)`). Instead, on the full cusp space `M`, intersect
`⋂_ℓ ker(T_𝔭 − c_ℓ) mod p` over **inert** primes `ℓ`, where the induced
structure gives the single linear condition
`c_ℓ = a_𝔭(f) = −[T² coefficient of the L-polynomial of A at ℓ] (mod p)`.
(At inert primes this value lies in `F_p`, even though split-prime eigenvalues
lie in `F_{p²}`.) Building `M` is fast (≤ 5 s even at dim 17550); the cost is the
`F_p` Hecke matrices and kernels.

---

## Feasibility wall (why 50 are out of reach)

Space dimension grows ~linearly in the level norm, with a field-dependent slope:
`Q(√2)` ≈ `0.072 × norm`, `Q(√3)` ≈ `0.173 × norm`. Dense `F_p` linear algebra is
feasible to dim ~30000 (≈ 1 GB per Hecke matrix). This cuts a clean line:

- **reachable:** mod-5 up to norm ~405000 over `Q(√2)` / ~101250 over `Q(√3)`;
  the three small-to-mid mod-7 entries.
- **out of reach:** the `serre = 2304`, `20736`, … buckets (real norms
  1.4M – 13M, dims up to 2.25M), plus the two large mod-7 entries (idx 73 norm
  40500 already dim 113398; idx 74 norm 162000). No mathematical shortcut removes
  the ×625, so these need a **custom sparse mod-p Hecke implementation**
  (Dembélé–Voight over `F_p` with the inert-kernel conditions) to go further.

---

## Mathematical takeaways relevant to Goal 2 (proving modularity)

- The **mod-5** forms have **small Hecke fields (degree 2–4)** ⟹ their attached
  `B_f` is a low-dimensional RM abelian variety over `F` (an abelian surface for
  degree 2). These are the realistic candidates for *explicitly constructing* `B`
  and comparing `A[5]` with `B[5]`.
- The **mod-7 idx-33** form has a **degree-18** Hecke field ⟹ `B_f` is
  18-dimensional, **not** an elliptic curve; the "`E_f`, `B = Res(E_f)`" plan
  does not apply there. Compare `A[7]` against `Ind ρ̄_f` directly instead.
- In all cases `A[p]` is **induced**, so its mod-`p` image is *imprimitive*
  (inside the normalizer of `GL₂×GL₂` in `GSp₄`). This is exactly the regime
  where the standard "big/enormous residual image" hypotheses of BCGP
  automorphy-lifting need care — **worth confirming which lifting statement
  covers the induced case** before claiming modularity of `A`.

## Artifacts

- `kernel_sweep.m` → `kernel_sweep_out.txt` — the control-validated sweep (all 25).
- `goal1_sweep.m` → `goal1_results.txt` — newform decomposition (Hecke degrees, small levels).
- `weilres_qsqrt10_*.m` — the idx-33 deep dive (trace + full char-poly confirmation).
- `goal1_data.m` — `sorted_output.jsonl` parsed to Magma (sorted by level).
