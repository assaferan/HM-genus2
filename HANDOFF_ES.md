# Handoff — analytic construction of B_f via EichlerShimuraHMF (2026-07-17)

Continues the Goal-2 (mod-5 idx-5) thread: get an **explicit genus-2 curve for `B_f`**, the
abelian surface over `F=Q(√2)` with RM by `Q(√3)` attached to the Hilbert newform
`f = 2.2.8.1-5000.1-j` (orbit 6, level norm 5000, weight [2,2], Hecke field `x²+2x−2`).

This session **closed the tooling gap** and computed `B_f`'s period matrix. What remains is a
precision (eigenvalue-count) grind, now running in the background and pushed to a branch.

---

## TL;DR status

1. **`B_f` is geometrically simple** — `End_Q̄⁰(B_f)=Q(√3)`, a *generic* point of the Humbert
   surface H₁₂ (NOT a badlocus building block — earlier memory claim corrected). Verified in
   `geom_endo_Bf.m`. ⇒ the point IS reachable; the obstruction is quantitative, not structural.
2. **Arithmetic H₁₂ CRT route: exhausted.** Height of the moduli point `> 2.57×10⁷`; brute-force
   CRT caps ~2.5×10⁷. See `geom_endo_Bf_notes.md` and `goal2_mod5_solve2.m`.
3. **Analytic route WORKS.** Using `edgarcosta/EichlerShimuraHMF` (Oda periods via twisted
   L-values) we computed `B_f`'s small period matrix:

   ```
   tau ≈ [ 1.16199 i   1.58643 i ]     Im(tau) positive definite (det ≈ 1.53):
         [ 1.58643 i   3.48608 i ]     a valid point of the Siegel upper half-space H_2
   ```
   at **~4-digit precision** (limited by the eigenvalue bound norm 10103). Script:
   `goal2_mod5_es_tau.m` → `goal2_mod5_es_tau_out.txt`.
4. **Remaining work = precision.** Recognizing the Igusa invariants over `Q(√2)` (→ the curve)
   needs ~30–60 digits ⇒ Hecke eigenvalues of `f` to norm ~10⁵. Expensive for level 5000
   (~0.8 s/prime, growing with `Norm(P)`), parallelizable, overnight-to-multiday. **Running now
   in the background, pushing to branch `es-highprec`** (see below).

---

## How to pick up on another machine

```bash
git fetch origin
git checkout es-highprec            # scripts, notes, and the growing eigenvalue cache
cat HANDOFF_ES.md
bash es_setup.sh                    # clone ../EichlerShimuraHMF + ../CHIMP siblings, register the form
tail -f es_eig/PROGRESS.md          # how far the eigenvalue precompute has gotten
```

Dependencies (siblings of this repo): **`../CHIMP`** (checked out, its spec attached) and
**`../EichlerShimuraHMF`** (cloned by `es_setup.sh`, which also re-applies the one-line
registration of our form in its `src/Labels.m`).

Once eigenvalues reach a higher norm `N` (check `es_eig/PROGRESS.md`), recompute the period
matrix at higher precision:

```bash
magma B:=60 MAXN:=<N> goal2_mod5_es_tau.m      # tau to ~ N/√(1.57e7)·1.5 digits
```

To get the **curve**: raise precision until `tau` is stable to ~40 digits, then feed it to the
Igusa/reconstruction (this repo's `Genus2Curve.m: IgusaInvariantsInK` / CHIMP `ReconstructCurve`)
and recognize the invariants over `Q(√2)`. VERIFY by matching `B_f`'s L-poly fingerprint (the
`<s,n>` table in `goal2_mod5_Bf_out.txt`).

---

## The five fixes that made ES work (do not re-discover)

All in `goal2_mod5_es_tau.m` (and `goal2_mod5_es_period.m`):

1. **In-process L-values.** The package's `ComputeSpecialValues` uses `ParallelPipe`; its spawned
   subprocesses fail here (`eval ... bad syntax`). Call `LSeriesTwisted` / `LMFDBTwistedLvalue`
   **directly** in one process (each L-value ~10 s; FE error ~1e-16, correct).
2. **`KnownConductor`** passed to `LSeriesTwisted` = `disc(F)²·Norm(level)·Norm(cond χ)²`
   = `320000·Norm(cond χ)²`. This **bypasses `GuessConductor`**, which needs far more eigenvalues
   than the L-value itself (it was the "Not enough eigenvalues to pin down the conductor" error).
3. **Coprime-conductor characters only** (`gcd(Norm cond χ, 10)=1`), so the `KnownConductor`
   formula is valid (χ ramified at 2 or 5 needs a different conductor).
4. **Polarization ideals `A=(1), B=(√3)`** (the ramified prime above 3). `Q(√3)`'s different
   `(2√3)` is mixed-sign, so the 17T7 example's `(1,1)` fails "`A·B·Different` narrowly principal".
5. **`tau` from three components** `Ω_{++}, Ω_{+-}, Ω_{-+}` — moduli point
   `z0 = [Ω_{-+}/Ω_{++} per Hecke embedding]`, then `SmallPeriodMatrix(z0, A, B)`. Skip
   `FixModuliPoint` (z0 already has positive imaginary parts).

### Known limitation to resolve for the CORRECT curve
`Ω_{--}` (sign `[-1,-1]`) currently errors: its smallest coprime character has conductor 9, whose
twisted L-conductor `320000·81 ≈ 2.6×10⁷` gives `MaxPrecision ≈ 1` digit at maxn=10103, and
`Evaluate` fails at that precision. `Ω_{--}` is NOT needed for `tau` itself, BUT it IS needed for
the correct *scaling* of the moduli point (the `cross_prod = Ω_{++}Ω_{--}/(Ω_{+-}Ω_{-+})` in
`PossibleModuliPoints`, which fixes a unit/rational rescaling → the principal polarization → the
right curve). Do NOT fabricate `Ω_{--}` from the relation — that makes `cross_prod ≡ −1` and
loses the moduli information. Fix = enough eigenvalues that the cond-9 `[-1,-1]` twist evaluates
to good precision (maxn ≳ 20000 for it to evaluate at all; more for precision).

**Precision rule of thumb:** `digits ≈ 1.5 · maxn / √(twist-conductor)`. Worst twist here has
conductor `320000·81 ≈ 2.6×10⁷` (√ ≈ 5100), so ~40 digits ⇒ maxn ≈ 1.4×10⁵.

---

## Eigenvalue precompute (the long pole)

- `f`'s Hecke eigenvalues cache to `es_eig/2.2.8.1-5000.1-j.txt` (LMFDB label : eltseq lines).
  Loaded into Magma's internal cache by `LMFDBHMFwithEigenvalues` ⇒ subsequent `HeckeEigenvalue`
  calls are instant. Currently **contiguous to norm 10103** (1230 primes).
- Cost: ~0.8 s/prime for small primes, growing with `Norm(P)` (T_P neighbor count). Reaching
  norm 10⁵ is overnight-to-multiday even parallel.
- **Generation:**
  - `es_eigenvalues.m` — single-process, RESUMABLE (loads existing, computes only missing up to
    `MAXN`), checkpoints every 250 primes. `magma MAXN:=30000 es_eigenvalues.m`.
  - `run_es_eig.sh MAXN NCHUNK` — parallel, interleaved, checkpointed chunks + merge (resumable
    via the base file). NOTE: 12 chunks thrash on redundant `NewformDecomposition`; use ≤ 8.
  - `hp_run.sh` — the background driver now running: loops over increasing norm targets, merges,
    and **commits + pushes `es_eig/*.txt` to branch `es-highprec`** after each milestone. Log:
    `hp_run.log`; human status: `es_eig/PROGRESS.md`.

---

## Key files (this session)

| File | Role |
|---|---|
| `geom_endo_Bf.m` | proves `B_f` geometrically simple (no self-inner-twist) |
| `geom_endo_Bf_notes.md` | full write-up: simple ⇒ generic H₁₂; height bound; g-route wall; ES success |
| `goal2_mod5_es_tau.m` | **THE working period-matrix driver** (all 5 fixes) |
| `goal2_mod5_es_period.m` | fuller driver (also tries the 4th sign class; verbose) |
| `es_register.m` | computes the registration data for `Labels.m` (level gens, trace fingerprint) |
| `es_eigenvalues.m`, `es_eig_chunk.m`, `run_es_eig.sh` | eigenvalue generation |
| `es_setup.sh` | reproduce the ES setup on a fresh machine |
| `hp_run.sh` | background high-precision loop + push to `es-highprec` |
| `es_eig/2.2.8.1-5000.1-j.txt` | the eigenvalue cache (grows) |
| `goal2_mod5_solve*.m`, `mcache.m` | the (exhausted) arithmetic H₁₂ CRT route |

Arithmetic/g-route/earlier Goal-2 context: see `goal2_mod5_notes.md`, `geom_endo_Bf_notes.md`.

---

## Standing result (unchanged, solid)

GRH-conditional modularity of the idx-5 (and idx-33) abelian surfaces is proven. The explicit
equation for `B_f` is now **reachable** (period matrix computed; only precision remains) — the
construction problem is no longer conceptually or tooling-blocked.
