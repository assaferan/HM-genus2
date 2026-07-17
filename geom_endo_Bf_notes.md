# `B_f` is geometrically simple — a generic H₁₂ point (2026-07-16)

**Context.** For the mod-5 idx-5 example, `B_f` is the abelian surface over `F=Q(√2)`
with RM by `K_f=Q(√3)` attached to the Hilbert newform `f` (orbit 6, level norm 5000,
weight `[2,2]`), LMFDB `2.2.8.1-5000.1-j`. The construction problem is to write down an
explicit genus-2 curve for it.

**What was previously concluded (and is now corrected).** The earlier H₁₂ reconstruction
(`goal2_mod5_loo2.m`) failed (best leave-one-out score 3/7 = chance), and this was read as
a *structural* obstruction: "the inner twist `a_{Pσ}=τ(a_P)` ⟹ (Ribet) `End_Q̄(B_f) ⊋
Z[√3]` ⟹ `B_f` sits on the H₁₂ extra-endomorphism badlocus ⟹ unreachable by the generic
2-parameter model." That chain conflates two different things.

**The distinction.**
- `a_{Pσ}=τ(a_P)` with `σ ∈ Gal(F/Q)` is the **Gal(F/Q)-conjugation** relation: `f^σ=f^τ`.
  It says `B_f` is a **Q-building block** (its Galois conjugate is isogenous) — a statement
  about the field of moduli / descent to `Q`. It does **not** create geometric endomorphisms.
- Extra **geometric** endomorphisms require a **self-inner-twist over `F`**: a *finite-order*
  Hecke character `χ` of `F` with `a_P·χ(P)=τ(a_P)`, i.e. `χ(P)=τ(a_P)/a_P` a root of unity
  at every good `P`.

**Computation (`geom_endo_Bf.m` → `geom_endo_Bf_out.txt`).** Over good primes to norm 80:
- `τ(a_P)/a_P` is an **infinite-order unit** at the split primes (`-2+√3` at 7, `-4√3-11`
  at 17, `(7√3+33)/23` at 23, …) — **never uniformly a root of unity**.
- `f` is **non-CM**: `a_P=0` at only 1 of 29 good primes.

So there is **no finite-order self-inner-twist**, hence

> **`End_Q̄⁰(B_f) = Q(√3)` — `B_f` is geometrically simple, a *generic* point of the
> Humbert surface H₁₂, not a badlocus point.**

Consistency check: `Res_{F/Q}(B_f) ~ A_g` is a simple degree-4-Hecke fourfold, and over
`Q̄`, `A_g ~ B_f²` with `End_Q̄⁰(A_g)=M₂(Q(√3))` (dim 8) — exactly what a geometrically
simple RM-by-√3 surface gives.

**Consequence for the construction.** The `loo2` failure was **not** a structural exclusion
from the `M_ℓ`. Since `B_f` reduces to a smooth generic H₁₂ point at every good `ℓ`, its
residue **is** in every `M_ℓ`, so reconstruction is possible *in principle*. The failure is a
**height wall**: `height(e,f)` exceeds the ~`10⁵` rational-reconstruction reach, so
rational-reconstruct returned a wrong small-height point (matches the CRT primes, fails fresh
ones → the chance-level 3/7).

**Net.** The H₁₂ **arithmetic** route is height-limited, not dead — pushing to higher height
with more/larger split primes and a better combinatorial `M_ℓ` solver can reach the true
point. The **analytic** route (period matrix of `B_f` → Igusa → `ReconstructCurve`) is the
height-independent alternative, still blocked only by the absence of a Hilbert-modular-form
period routine (or, equivalently, by needing the classical descended form `g` at the heavy
level 320000).

Scripts: `geom_endo_Bf.m`, `geom_endo_Bf_out.txt`.

## Arithmetic H₁₂ reconstruction, pushed (2026-07-16)

Acting on "generic H₁₂ point ⟹ reachable," I rebuilt the CRT reconstruction:

- **Fast `M_ℓ` builder**, cached for 20 split primes ≤199 (`goal2_mod5_mbuild.m` →
  `mcache.m`). `|M_ℓ|` ranges 2–172; the cheap ones (small `|M|`, good for CRT) are
  `7(2), 23(2), 17(6), 31(10), 41(10), 113(14), 137(16), 193(16), 73(26), 127(26)`.
- **Knapsack solver** (`goal2_mod5_solve.m`, `goal2_mod5_solve2.m`): choose a prime subset
  maximizing modulus `∏ℓ` (⟹ rational-reconstruction reach `~√(M/2)`) under a brute-force
  budget `∏|M_ℓ|`; verify candidates by membership at the remaining primes, allowing a few
  bad-reduction primes; plus **leave-one-out** over the CRT primes (cheap: `∑1/|M_ℓ|`).
- **Positive control** (`goal2_mod5_control.m`): a known `(e,f)=(3,-2)` has its true residue
  in `M_ℓ` at every good prime and reconstructs — **the builder is validated**, so `B_f`
  null results are genuine reach limits, not bugs.

**Results.**
| run | reach | combos | best score |
|---|---|---|---|
| 8 primes ≤193 | 2.28×10⁶ | 8.6M | 4/12 (chance) |
| + leave-one-out | ≤2.28×10⁶ | ~26M | 4 (chance) |
| 9 primes (+127) | **2.57×10⁷** | 224M (62 min) | 4/11 (chance) |

⟹ **`height(e,f) > 2.57×10⁷`.** The point is genuinely high.

**Ceiling.** Reach `~√(∏ℓ)` but cost `~∏|M_ℓ|`, so brute-force CRT caps near ~2.5×10⁷ with
primes ≤199. More cheap primes are impractical (`M_ℓ` build is ~1000 s/prime at ℓ~230 via the
F_ℓ² scan, and small-`|M|` large primes are rare: 223, 233 gave `|M|`=30, 372). The tempting
"membership-filter instead of rational reconstruction to dodge the height bound" fails — the
*integer* CRT residue reduces correctly only at the CRT primes, so `RationalReconstruction`
(needing `modulus > 2·height²`) is unavoidable.

**Bottom line.** `B_f` is a generic, reachable H₁₂ point, but its moduli coordinates have
height `> 2.6×10⁷`; blind arithmetic CRT is exhausted. Finishing needs a **real-number seed**
for `(e,f)` (even ~10–20 digits) to switch from blind combinatorics to direct LLL/recognition
— i.e. the analytic period input (Hilbert-form periods, tooling-blocked; or the classical `g`
at level 320000 via modular symbols, heavy but standard).

Scripts: `goal2_mod5_mbuild.m`, `goal2_mod5_solve.m`, `goal2_mod5_solve2.m`,
`goal2_mod5_control.m`, `goal2_mod5_extend.m`, `mcache.m`.

## Analytic route WORKS: B_f period matrix via EichlerShimuraHMF (2026-07-17)

The tooling gap is closed. `edgarcosta/EichlerShimuraHMF` (Oda periods of Hilbert modular
abelian varieties via twisted L-values; dep = CHIMP only) handles our exact case. Registered
`2.2.8.1-5000.1-j` in its `Labels.m` (dim 2, Hecke field `x²+2x−2`, level `50√2`), precomputed
Hecke eigenvalues (contiguous to norm 10103), and reconstructed **`B_f`'s small period matrix**:

    tau ≈ [ 1.16199 i   1.58643 i ]      (Im tau positive definite, det ≈ 1.53:
          [ 1.58643 i   3.48608 i ]       a valid point of the Siegel upper half-space H_2)

Fixes required (all in `goal2_mod5_es_tau.m`):
1. **In-process L-values** — the package's `ComputeSpecialValues` uses `ParallelPipe`, whose
   spawned subprocesses fail here; call `LMFDBTwistedLvalue`/`LSeriesTwisted` directly instead.
2. **`KnownConductor`** = `disc(F)²·Norm(level)·Norm(cond χ)²` (= `320000·Norm(cond χ)²`) passed to
   `LSeriesTwisted` — bypasses `GuessConductor`, which needs far more eigenvalues than the value.
3. **Coprime-conductor characters only** (so the `KnownConductor` formula is valid).
4. **Polarization ideals** `A=(1), B=(√3)`: `Q(√3)`'s different `(2√3)` is mixed-sign, so the
   example's `(1,1)` fails the "`A·B·Different` narrowly principal" requirement.
5. `tau` comes from the **three** components `Ω_{++}, Ω_{+-}, Ω_{-+}` (moduli point
   `z0 = Ω_{-+}/Ω_{++}` per Hecke embedding). `Ω_{--}` errors (its doubly-ramified twist has too
   little precision at maxn=10103) and only fixes a unit scaling, so it is not needed for `tau`.

**Precision is ~4 digits**, limited by the cond-7 twisted L-values at maxn=10103 (precision ~
maxn/√(twist-conductor)). Recognizing the genus-2 curve over `Q(√2)` needs ~30–60 digits ⇒
eigenvalues to norm ~10⁵, an expensive precompute for our level-5000 form (~0.8 s/prime growing
with `Norm(P)`; parallelizable but overnight+). Scripts: `goal2_mod5_es_tau.m`,
`goal2_mod5_es_period.m`, `es_eigenvalues.m`, `es_register.m`, `es_eig/` (eigenvalue cache).

## g-seed route scoped (2026-07-16): Magma modular symbols blocked at dim 48004

To get a real-number seed for `(e,f)` via the classical descended form `g` (weight 2, level
`N=320000=2⁹·5⁴`, nebentypus `χ₈` = quadratic character of `Q(√2)`, degree-4 Hecke field
`K_g=Q(√2,√3)`, with `A_g=Res_{Q(√2)/Q}(B_f)`):

- `dim S_2(N,χ₈)` FULL = **47520**, NEW = **7680** (formula, instant).
- `ModularSymbols(χ₈, 2, +1)` ambient **builds in 93 s** (dim 48004) — fine.
- **But `CuspidalSubspace` at dim 48004 does not finish in 15 min** and is killed. Cutting
  `g` (needs the cuspidal/new subspace or a char-0 kernel) and then computing periods are
  both heavier still. This is the same dim ~30k–60k char-0-linear-algebra wall hit in Goal-1.

⟹ **The Magma modular-symbols `g`-route is not viable** at this level. The bypass is a
from-scratch numerical Eichler integral using `a_n(g)` from `B_f`'s Euler product (which we
can evaluate at any prime from `f`), folded via the functional equation (conductor 320000) to
~15-digit precision — enough to seed LLL against the height-`2.6×10⁷` rational. That is a
substantial implementation (get `a_n(g)∈K_g` from the automorphic induction, numerical
periods, then split off `B_f`'s 2×2 `τ` via the two `√2→+` conjugate embeddings). Frontier
work; not a quick win. Scripts: `goal2_mod5_g_scope.m`.
