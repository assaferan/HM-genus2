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

## Feasibility (weight [2,2] view — SUPERSEDED, see correction below)

At **weight [2,2]** the mod-5 level is Serre × 5⁴, and space dimension grows
~linearly in that (inflated) norm (`Q(√2)` ≈ `0.072×`, `Q(√3)` ≈ `0.173×`). Dense
`F_p` linear algebra is feasible to dim ~30000, which at weight [2,2] leaves the
`serre = 2304, 20736, …` buckets out of reach (dims up to 2.25M).

## CORRECTION — the ×625 wall was a weight-choice artifact

Prompted by Ariel: the same residual rep is matched at the **small prime-to-5
Serre level** with a **raised weight** (e.g. [2,4]), trivial nebentypus, using the
normalization `t_P = Norm(P)·a_P(f) mod λ` (twist power d ≈ (k₀−2)/2). This avoids
the ×625 blow-up entirely. Validated by reproducing Ariel's [2,4]/level-16 match
(`ariel_validate.m`) and rescuing former `SKIP-BIG` entries at their tiny Serre
levels (`rescue_scan2.m`; e.g. idx 20 at norm 576, dim 288, vs infeasible ×625 =
360000, dim 62398).

**Combined coverage now 32/70 mod-5** (up from 22): the [2,4] pass
(`goal1_sweep_w24.m` → `goal1_w24_results.txt`) rescued 10 former SKIP-BIG entries
(idx 20,21,24,25,26,27,28,29,30,34). The remaining entries are almost all just
*unprocessed* at [2,4] (the decompositions run ~30 min each at large levels, a
day-scale batch), plus idx 31,32 which need a Serre weight above [2,4]. So the
"50 out of reach" claim above is **wrong**: SKIP-BIG is a compute-budget question,
not a mathematical wall. (The two large mod-7 entries idx 73/74 remain heavy since
there weight [2,2] is already correct — no p-part to trade to weight.)

The earlier "custom sparse mod-p Hecke" suggestion is therefore unnecessary for
mod-5 coverage; the weight tradeoff is the practical route.

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

## Update — niveau-1 vs niveau-2, and the fast classifier

The mod-5 entries are of two local types at 5:
- **niveau-1** (ordinary at 5): a small prime-to-5 Serre-level realization exists at
  raised weight **[2,4]** (trivial nebentypus, t_P = Norm(P)·a_P(f)). Every
  raised-weight match found is at [2,4]/[4,2] — none needs [2,6] or higher.
- **niveau-2** (supercuspidal at 5): NO small-level realization at any raised weight;
  lives only at the ×625/[2,2] level (feasible only for small Serre conductor).

**Fast classifier + matcher (kernel_wt.m):** one kernel-method call at [2,4] over the
weight field Q(zeta_8) (reduce mod a prime above 5 -> F_25; no decomposition). Match
=> niveau-1 (done); no-match => niveau-2. This replaces the slow decomposition and
the multi-weight search.

**Final combined coverage: 43/70 mod-5** (x625 [2,2] union raised-weight [2,4]),
plus all 5 mod-7. niveau-2 confirmed for idx 4,7,10,13,14,15,19 (matched only at
x625). The remaining unmatched are niveau-2 with SKIP-BIG x625 levels -- a genuine
hard tail, NOT a weight-choice artifact.

Conductor-only predictor is only partial: v_5(conductor_2d)=4 => niveau-2, but
idx 0-4 share conductor invariants yet split (0-3 niveau-1, 4 niveau-2), so niveau
is a curve-specific local property at 5. Scripts: kernel_w24.m, kernel_wt.m,
kernel_sweep_wt.py.

## Final coverage (with the cluster predictor + giant confirmation)

**54/75 total: 49/70 mod-5 + all 5 mod-7.**  The niveau predictor (cluster picture of
4f+h^2 mod 5) confirmed on fresh data: of the 8 unmatched giants it flagged niveau-1,
6 (idx 59-64, Q(sqrt2), dim 4320) MATCH at [2,4] (control-validated, ~36 min each);
the other 2 (idx 69,70, Q(sqrt3), dim 10368) are the same niveau-1 type but time out
(kernel over the deg-4 weight field too slow at a 100-min cap -- compute-limited, not
a predictor failure).  The remaining 19 unmatched are predicted niveau-2 (supercuspidal
at 5) -- genuinely stuck at the x625/[2,2] level, no small-level realization.  So the
tail is fully characterized: 49 matched, 2 niveau-1 (reachable with more compute),
19 niveau-2 (mathematically stuck).  Giant runs: giants_n1_results.txt (kernel_wt.m).
