# Handoff — mod-11 Serre modularity evidence over Q(√5)

Orientation for another Claude instance picking up the **Serre-modularity sub-project** on
branch `serre-mod11-evidence` (PR #2). For the broader HM-genus2 project see `HANDOFF.md`,
`CLAUDE.md`, `README.md`. Full technical write-up: **`results_mod11_serre.md`** and the LaTeX
note **`serre_mod11_note.tex`/`.pdf`**. All work committed + pushed.

## The goal

Exhibit computational evidence for **Serre's modularity conjecture over the real quadratic
field K = Q(√5)** (open there): take an abelian surface whose mod-11 Galois representation has
a 2-dimensional piece, and match that piece to a Hilbert modular form.

## The object: "Curve A"

Genus-2 curve `C/Q(√5)` from the `(1,11)` Gross–Popescu construction (a point on the `D₂`
sextic). Igusa invariants in `verify_curveA.m`/`results_mod11_serre.md`. Verified:
`End_{K̄}=Z` (geometrically simple, **no RM**), good reduction at both primes above 11,
multiplicative reduction at the prime `𝔭` of norm **66179**, and (model-)bad at the inert
prime `(2)`; reduction at 5 is the "model-degenerates" prime (uncertain).

## What is established (solid)

1. **`A[11]` splits `2+2` over K.** The Frobenius quartic mod 11 factors only as `[1,1,2]` or
   `[1,1,1,1]` over 290 primes (never `[2,2]/[1,3]/[4]`) — forced by the `(1,11)` polarization
   at 11 (11 = polarization degree), NOT by endomorphisms. So `A[11] = ρ₁ ⊕ ρ₂`, two
   irreducible 2-dim reps, `det = χ₁₁`, large image, not twists of each other, no 1-dim
   K-subquotient. (`factor_types.m`, `analyze_two_reps.m`, `check_endo.m`, `determine_lambda.m`.)
2. **Fingerprint.** `{tr ρ₁, tr ρ₂}` = roots of `y² + A y + (B−2N(P))` mod 11, where `h_P =
   x⁴+Ax³+Bx²+Cx+D`. The symmetric `(Tr,Nm)=(−A, B−2N)` is in **`mod11_trnm.txt`** (290 primes);
   raw L-polys in `mod11_frobdata.txt`. Built by `mod11_rep.m`.
3. **`ρ₁` is MODULAR.** Working entirely mod 11 (definite quaternion algebra; the Hilbert space
   is cheap — dim 1102 at norm 66179, builds instantly), the twist-robust kernel intersection
   `V = ⋂_P [ker(T_P²−Tr_P T_P+Nm_P) + ker(T_P²+Tr_P T_P+Nm_P)]` finds an eigensystem `σ₁` that
   is a **Hilbert newform of parallel weight 2 and level `𝔭` (norm 66179)** with
   `tr σ₁ = ε·tr ρ₁`. Verified `a_P² = tr(ρ₁)²` over **24 primes**. (`match_hmf_twist.m`,
   `extract_eigensystems.m`, `level_p_twist.m`.)
4. **The modular form gives the global labeling** of `ρ₁` vs `ρ₂` — no `K(L)`/Heisenberg
   computation needed. `tr ρ₁ := ` the root whose square is `a_P(σ₁)²`; `tr ρ₂` = the other root.
5. **The twist is identified:** `ε = χ_{K(√d)}`, `d = 66179·(5−√5)/2`, ramified at the prime
   above 5 and the primes above 66179, **unramified at 2**. Matches **17/18** non-degenerate
   primes. (`identify_twist2.m`, `verify_twist.m`, data in `eps_data.txt`.)

## Open items (the honest gaps)

- **`ρ₂`'s form has not been located.** Checked and ABSENT at every level tried: `𝔭`,
  `𝔭^τ` (the other prime above 66179), `𝔭·(2)`, and **`𝔭·(2)²` (dim 22059, run on lovelace —
  `dim V = 3` = `σ₁` oldforms only, no new system)**. Its traces are explicit (complementary
  roots) but its form is unknown. **Best next candidate: levels with ramification at 5** —
  `𝔭·(5)` (norm 330895), `𝔭·(2)·(5)`, `𝔭^τ·(5)` — consistent with the conductor subtlety below.
  These are dim ~10⁴–10⁵; use `lovelace` (see below). Reuse the `rho2_search.m` method
  (`ker(T_P−λ)` + restrict-to-cur).
- **One residual twist misfit:** `⟨59, rid 8⟩` (one prime above 59). It is GENUINE, not an
  extraction glitch (`σ₁` is a multiplicity-one newform at `𝔭`, so `a_P` is unambiguous). So
  either `ε` needs a slightly different conductor at that prime, or Curve A's L-poly at that
  prime has an error worth re-checking in `mod11_frobdata.txt`.
- **Conductor subtlety (unresolved):** `ε` looks ramified at `{5, 66179}` yet `σ₁` sits at
  level `𝔭` (no 5-part) — consistent only if `ρ₁` itself is ramified at 5 (plausible: 5 is the
  bad-model prime). The exact `σ₁ = ρ₁⊗ε` conductor bookkeeping is not nailed. Worth: determine
  Curve A's actual reduction/conductor at 5 and at 2 (Magma's genus-2 conductor fails at these
  over number fields — may need a regular-model or cluster-picture approach).

## Key files (this sub-project)

| File | Role |
|---|---|
| `results_mod11_serre.md` | full technical write-up (READ FIRST) |
| `serre_mod11_note.tex/.pdf` | 4-page math note (honest, computational-evidence framing) |
| `mod11_rep.m` | reduce Curve A mod P, build `mod11_frobdata.txt`, `mod11_apdata.txt` |
| `mod11_trnm.txt` | the `(Tr,Nm)` fingerprint of `{ρ₁,ρ₂}` — the matching target |
| `match_hmf_twist.m` | twist-robust mod-11 match (the dim-2 → `σ₁` result) |
| `level_p_twist.m` | shows `σ₁` is a newform at level `𝔭`; tries `𝔭^τ` |
| `extract_eigensystems.m` | extract `σ`'s `a_P`, compute `ε(P)`, save `eps_data.txt` |
| `identify_twist2.m`, `verify_twist.m` | identify `ε = χ_{K(√d)}` |
| `rho2_search.m` | the `ρ₂` search at `𝔭·(2)²` (dim 22059) — built for lovelace |
| `factor_types.m`, `analyze_two_reps.m`, `check_endo.m` | structural findings |

Earlier matcher iterations (`hmf_match.m`, `match_hmf_mod11.m`, `build_match_hmf.m`,
`debug_*.m`, `determine_lambda.m`) are kept for provenance. `*_progress.txt` are gitignored.

## Compute environment

- **Local mac:** Magma **V2.29-3**. CHIMP at `../CHIMP`. If CHIMP fails to attach with a
  "compiled for newer version" error, clear caches: `find ../CHIMP -name '*.sig' -delete`.
- **Remote `lovelace` (use for heavy jobs):** see the **`lovelace-remote-compute`** memory.
  256 cores, ~2 TB RAM, `ssh lovelace` (key auth). CHIMP at `~/GitHub/CHIMP`; clone the repo
  as a sibling: `cd ~/GitHub && git clone https://github.com/assaferan/HM-genus2.git`. Launch
  detached: `ssh lovelace 'cd ~/GitHub/HM-genus2 && setsid bash -c "magma JOB.m > JOB.log 2>&1"
  < /dev/null &'`; have the script `PrintFile` to a `*_progress.txt` and poll it via ssh.
  GOTCHA: `pkill -f "magma JOB"` to stop a run; two runs writing the same progress file
  interleave (start with `Overwrite:=true` and use a fresh name per run).

## Suggested next steps (in order)

1. **Hunt `ρ₂` at 5-ramified levels** on lovelace: `𝔭·(5)`, `𝔭·(2)·(5)`, `𝔭^τ·(5)` with the
   `rho2_search.m` method. A new eigensystem with `a_P` = the complementary roots would exhibit
   `ρ₂` modular and complete the picture (both pieces modular ⇒ `A[11]` fully modular).
2. **Resolve the conductor story:** compute Curve A's reduction at 5 and 2 (regular model /
   Liu / cluster pictures), confirm `ρ₁` ramified at 5, and reconcile `σ₁ ∈ S_2(𝔭)` with
   `ε` ramified at `{5,66179}`.
3. **Clear / explain `⟨59,8⟩`:** re-verify Curve A's L-polynomial there in `mod11_frobdata.txt`.
4. **Extend verification** of `σ₁ ↔ ρ₁` to ~40+ primes for a tighter statement.
