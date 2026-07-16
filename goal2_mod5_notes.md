# Goal-2 for a mod-5 example (idx 5) — the tractable analog of idx 33

**Motivation.** The mod-7 flagship (idx 33) gave a GRH-conditional modularity
certificate but its unconditional route was blocked because the matching newform
has a **degree-18 Hecke field**, so `B_f` is an irreducible dim-18 abelian
variety — not constructible. The mod-5 matches all have **degree-2 Hecke fields**,
so `B_f` is an **abelian surface** and `Res_{F/Q}(B_f)` is 4-dimensional, exactly
the regime of van Bommel–Costa–Elkies–Keller–Schiavone–Voight (arXiv:2411.07857).
This note runs the idx-33 Goal-2 pipeline on a mod-5 example to see how far it goes.

**Example (idx 5).** `F = Q(√2)`, `p = 5`,
`C: y² = 5 − 10x + 10x² + 5x⁴ + 2x⁵`, `cond(A) = 1 600 000 = 2⁹·5⁵`.
Goal-1 match: **orbit 6** of `HilbertCuspForms(F, 𝔫, [2,2])`,
`𝔫 = P₂³·P₅²` (norm 5000, P₂ the ramified prime over 2, P₅ the inert prime over 5),
Hecke field degree 2, reduced at `λ|5` with residue field **F₂₅**.
Note 5 is a **bad** prime here (5 | cond), so the level carries `P₅²` — unlike
idx 33 where 7 was good and the level was the Serre ideal exactly.

## Step 1 — image of A[5] (`goal2_mod5_idx5.m`, DONE)

From A's char polynomial of Frobenius mod 5, primes up to 2000:
- **Induced**: 155/155 primes inert in F have `a_ℓ ≡ 0 (mod 5)` (the swap-coset
  signature). 0 failures.
- **Irreducible**: 121/146 split primes give the σ-trace pair as a conjugate pair
  in `F₂₅∖F₅` (nonsplit-semisimple), so `σ ≇ σᵗ` ⟹ `Ind σ` irreducible.
- **Big image**: an order-5 transvection is present and there are 12 distinct
  projective traces `t²/ℓ`; by **Dickson**, `im(σ) ⊇ SL₂(F₂₅)`. (The raw
  "tr = 0 at 60.8%" line in the log is dominated by the 155 forced-zero *inert*
  primes; among split primes only 28/146 ≈ 19% — the transvection alone excludes
  the dihedral/normalizer-of-torus case.)

## Step 4 — GRH-conditional certificate (`goal2_mod5_idx5.m` + `..._cert2.m`, DONE)

`A[5]` and `Ind(ρ̄_f)` are both 4-dim, irreducible, unramified outside `S = {2,5}`,
with Artin conductor supported there. By Brauer–Nesbitt they are isomorphic iff
their Frobenius traces agree; under GRH the conductor-based effective bound on the
least distinguishing prime is `O((log cond)²)` — here `(log 1.6e6)² ≈ 205`, so a
pass to `10³–10⁴` suffices, **independent** of the huge splitting field.

**Verified** `tr A[5](Frob_ℓ) = tr Ind(ρ̄_f)(Frob_ℓ)` for **all good primes ℓ < 12000**
(≈1437 primes; A-side `a_ℓ(A) mod 5`, Ind-side `a_𝔭(f)+a_𝔭'(f)` at split ℓ and 0 at
inert ℓ), **0 disagreements** (run 1: ℓ<5003, 668 primes; run 2: 5000≤ℓ<12000, 769
primes). This exceeds idx-33's 1225-primes-to-10⁴ pass.

**⇒ Under GRH, `A[5] ≅ Ind(ρ̄_f)`**, so `A[5]` is residually modular. With GL₂/F
automorphy lifting (as in idx 33, Ariel's route) this yields GRH-conditional
modularity of A.

## Where mod 5 genuinely helps — and where it does not

The torsion-field route to drop GRH (van Bommel et al.) has **two walls**; mod 5
removes one:

- **B_f side — REMOVED.** `B_f` is a degree-2-RM **abelian surface** over `Q(√2)`,
  `Res_{F/Q}(B_f)` is dim 4 — constructible in principle and within the paper's
  precedent. The *modular object is now concretely exhibitable* (idx 33's dim-18
  `B_f` was not).
- **A side — STILL PRESENT.** The projective 5-torsion resolvent `F_A` blows up
  because **End(A) = Z** (no RM to expose the F₂₅-rational structure) — a property
  of A, independent of reducing mod 5 vs mod 7. Same wall as idx 33.

Also BCGP Prop 10.1.3's clean *p=5* induced case wants **residue degree 1** (F₅);
our form has 5 inert in its Hecke field `Q(√2)` ⟹ residue field **F₂₅**, the same
residue-degree-2 subtlety as idx-33's F₄₉.

**Net.** mod 5 buys a cleaner, faster certificate and — the real new thing — a
*constructible* modular abelian surface `B_f`. It does **not** by itself deliver
unconditional modularity: the A-side End=Z wall and the residue-degree-2 subtlety
remain. Next step: actually construct `B_f`.

## Constructing B_f (`goal2_mod5_Bf.m`, DONE up to identification)

**B_f identified as a named modular object.** The orbit-6 eigenform is **LMFDB
Hilbert newform `2.2.8.1-5000.1-j`** (verified: Hecke field `x²+2x−2 = Q(√3)`,
and `a_P` at norms 7/17/23 are `{−e−2, e}`, `{−2e+1, 2e+5}`, `{−e+6, e+8}`,
matching our Magma fingerprint exactly). So

> **B_f = the abelian surface over Q(√2) with RM by Q(√3) attached to
> `2.2.8.1-5000.1-j`**, conductor 𝔫 (norm 5000), GL₂-type.

(LMFDB form `-f` is `x²−2x−2` giving `1±√3` — a *different* newform; `-g,-h,-i`
have Hecke field Q(√2). Only `-j` matches A.)

**Exact Frobenius fingerprint** (`goal2_mod5_Bf_out.txt`): at a good prime P,
Frob on B_f has char poly `x² − a_P x + NP` with `a_P ∈ Z[√3]`, and the degree-4
Q-L-polynomial is `Norm_{Q(√3)/Q}(1 − a_P T + NP T²) = NP²T⁴ − NP·s T³ +
(2NP+n)T² − s T + 1`, `s = Tr(a_P)`, `n = Nm(a_P)`. Tabulated for 21 good primes;
at inert P (norm ℓ²) `a_P ∈ Z`, at split P the two primes give Galois-conjugate
`a_P` — the RM signature.

**No explicit model is tabulated anywhere** (LMFDB: "L-function not available",
no related genus-2 curve / abelian variety). Producing a defining equation needs
the **discriminant-12 Humbert/Hilbert-modular-surface parametrization** (Bruin–
Flynn–Shnidman `Y_(12)[√3]`; Elkies–Kumar `H_12`; generic models arXiv:2403.03191,
a 3-parameter family `y² = Nm_{L/K}φ(x)`, `L = k[r]/(r³−3(a²−3b²)r+2(a²−3b²))`,
with a conic obstruction) **plus a moduli-point search over Q(√2)** for the point
whose Jacobian reproduces the fingerprint above. This is the well-defined but heavy
next step; the explicit f₁₂ lives in that paper's electronic supplement.

**Why this is the payoff vs idx 33.** idx-33's B_f is an irreducible **dim-18**
variety — not in any parametrized family, not constructible. Here B_f is an
abelian surface lying in a **known 2-dimensional family** of explicitly
parametrized RM-by-√3 Jacobians, so an equation is in-principle reachable. That is
exactly the tractability gain the mod-5 route was expected to give.

## Chasing an explicit equation (`bfs_family.m`, `goal2_mod5_search*.m`)

**Family obtained and validated.** Transcribed the Bruin–Flynn–Shnidman explicit
family (arXiv:2102.04319, Section 3, "√3 over the ground field"):
`C_{a,b,c}: y² = G₁(x)² + λ₁H₁(x)³`, with `G₁,H₁,λ₁ ∈ Z[a,b,c]` (leading term of
G₁ is `b(a−c)²(a²−4ac−b²+c²)q₄²x³` — note the (a−c)², a transcription trap), q₄ =
Q(c,b,a), t₃ = Q(a,ζb,c)Q(a,ζ²b,c), Q(X,Y,Z)=X²+X(Y−Z)+(Y−Z)²−3XZ. Validated
against the paper's sample: (a,b,c)=(1,2,−1) reproduces `y²=8x⁵−3x⁴−2x³−7x²+4x+20`
(IsIsomorphic = true). The fast mod-ℓ point-counter is validated against Magma.
Machine-readable source: `math.huji.ac.il/~shnidman/BFScode.sage`.

**Moduli search — not yet pinned.** Target: (a:b:c) ∈ P²(Q(√2)) with Jac ≅ B_f,
filtered by |trace| (twist-invariant) at split primes 7,17,23,31,41,47 against the
fingerprint (targets s_ℓ = −2,6,14,6,−8,0). Searched: (i) blind box a,b,c∈Z[√2],
height ≤3 — 0 hits; (ii) {2,5}-S-unit set, c=1, a,b over 432 S-units — 0 hits;
(iii) blind box height ≤6 — running.

**The obstruction (why the clean shortcuts fail).** The BFS family requires full
**√3-level structure**. B_f arises from the mod-**5** representation, so it has no
reason to carry rational √3-level structure over Q(√2); hence B_f is a family
member only up to a **quadratic twist by some d ∈ Q(√2)\***, and d's conductor need
not be {2,5}. That breaks the S-unit shortcut (the untwisted C_{a,b,c} has bad
reduction at d's primes, so a,b,c are NOT {2,5}-units) and leaves the moduli
point's height unknown. |trace|-matching still finds it (twist-invariant), but only
once the box reaches the right height.

**DECISIVE DIAGNOSIS (`goal2_mod5_igusaset.m`): B_f is NOT a quadratic twist of any
BFS family member.** For each split prime ℓ we scanned ALL of P²(F_ℓ) and kept
curves whose Jacobian L-poly equals B_f's target L-poly *up to quadratic twist*
(L_target or its s→−s twin). Result:
  ℓ=17: 32 matching points (8 distinct Igusa tuples);
  ℓ=7, 23, 31: **0 matching points**.
If B_f were a quadratic twist of some C_{a,b,c}, its reduction mod every good prime
would be a matching family point — in particular a smooth F₇-point, of which there
are none. So no (a:b:c)∈P²(Q(√2)) has Jac(C_{a,b,c}) a quadratic twist of B_f.

**Root cause.** The BFS family requires full **√3-level structure** (the order-3
subgroups D₁,D₂ rational). B_f comes from the mod-**5** representation and carries
no such structure, and a quadratic twist cannot supply it. Relatedly, the modular
surface B_f attached to a Hilbert eigenform is **generally not principally
polarized** (BFS note their A/B_{a,b,c} are not PP), i.e. B_f need not be a
Jacobian at all — only isogenous to one. So the BFS (level-structured) family is
structurally the wrong parametrization here; the earlier height/S-unit search
failures were symptoms of this, not just large height.

**Correct tool + remaining work.** Use the **Humbert surface H₁₂** (Elkies–Kumar,
arXiv:1209.3527), which parametrizes RM-by-√3 Jacobians *without* level structure
via Igusa–Clebsch invariants (I₂:I₄:I₆:I₁₀) as rational functions of 2 parameters.
A PP genus-2 Jacobian isogenous to B_f (if one exists over Q(√2)) is a point of
H₁₂; pin it by matching L-polys / reconstructing its Igusa invariants (the ℓ=17
foothold gives 8 candidate Igusa tuples mod 17). Alternatively the **analytic
period** route: period matrix of B_f from the Hilbert eigenform → Igusa invariants
→ CHIMP `ReconstructCurve` (blocked only by the missing Hilbert-form-period
routine in Magma). Status: family + search infrastructure validated; the structural
obstruction to the BFS route is now understood; an explicit model needs H₁₂ or
analytic periods.

## The H₁₂ route (the correct family) — `goal2_mod5_h12.m`, `..._crtQ.m`, `..._curve.m`

**H₁₂ confirmed as the right family.** From the Elkies–Kumar ancillary data
(arXiv:1209.3527, `12/igusa12.txt`) the discriminant-12 Humbert surface has a clean
**2-parameter** model: with `A1=f−1, A=−(f⁴+15ef+9e)/3, B1=(2f³−2f²−3f+3e+3)/3,
B=(2f⁶−63ef³−81ef²−54e²)/27, B2=e³`, the Igusa–Clebsch invariants are
`[I2,I4,I6,I10] = [−24B1/A1, −12A, 96(A/A1)B1−36B, −4A1B2]` (bad locus
`e·f·(f²−1)·(e−f²)·(8f⁴−9ef²−3e)·(f⁶−f⁴−18ef²+27e²+16e)`). Unlike BFS, this needs
**no √3-level structure**, so it contains B_f's isogeny class. Per-prime scan
(`goal2_mod5_h12.m`): B_f's L-poly is matched (up to twist) at **every** split
prime 7,17,23,31 — the exact opposite of BFS's 0-matches. ✓

**Reconstruction machinery built and run** (`goal2_mod5_targets.m` extends the
fingerprint to 20 split primes; `goal2_mod5_mcount.m` finds small-|M| primes;
`goal2_mod5_crtQ.m` does the CRT). Key idea: at each split prime the true (e,f)
reduces into the small match-set M_ℓ; CRT + rational reconstruction across primes
recovers it. **Structural finding:** the two spurious survivors that first appeared
were *rational* — consistent with f having an **inner twist** (a_{P^σ}=τ(a_P) at
conjugate primes, seen in the fingerprint), i.e. B_f's field of moduli is Q, so
(e,f) ∈ Q² and the CRT needs only ∏|M_ℓ| combinations (not ∏|M_ℓ|²).

**Moduli point NOT yet reconstructed — the remaining wall.** Reconstruction was
pushed to height reach ~2730 (over Q(√2)) and **~10⁵ (rational, `crtQ` with primes
7,23,17,31,41,73,79)** — always **0 true survivors** (only twist/L-poly false
positives that die at fresh primes). A height >10⁵ for a conductor-5000 moduli
point is implausible, so the likely cause is a **badlocus collision**: B_f's
reduction lands on the H₁₂ extra-endomorphism / model-degeneration locus at one of
the small-|M| CRT primes (7,23,17,31), which *excludes* its residue from M_ℓ, so no
modulus can recover it. The clean fix is **leave-one-out** reconstruction (drop each
prime; the run omitting the bad ℓ₀ succeeds) — attempted but the large-prime F_ℓ²
match-set scans (|M| up to 116 at ℓ=97) are prohibitively slow; needs a faster
match-set builder (e.g. only over the ordinary/generic locus, or point-count L-poly
without full Mestre reconstruction). That is the concrete next step.

**RESOLUTION (`goal2_mod5_loo2.m`, the fast builder): B_f is a BUILDING BLOCK,
lying on the H₁₂ badlocus — so it is NOT reachable by the generic 2-parameter model.**
The fast match-set builder (trace-prune: cheap #C(F_ℓ) first, expensive #C(F_{ℓ²})
norm check only on the few trace-matches) cut match-set construction from ~700 s to
~19 s, making an exhaustive search feasible. Result:
 - rational CRT reconstruction over pool {7,23,17,31,41,73,79} reaches height 10⁵;
 - **drop-0, all drop-1, and drop-2** reconstructions were scored by how many of 7
   fresh primes each candidate matches — the **global best is 3/7 (chance level)**.
So B_f's moduli point is reconstructed by NO subset — not a height/single-badlocus
issue.

**Why:** the fingerprint shows an **inner twist** — `a_{Pσ} = τ(a_P)` (σ = Gal(Q(√2)/Q),
τ = Gal(Q(√3)/Q)) at *every* conjugate prime pair (e.g. −1∓√3 at the primes over 7).
By Ribet's theory of Q-surfaces/building blocks, a Hilbert eigenform with an inner
twist has `End_{Q̄}(B_f) ⊋ Z[√3]` (quaternionic or product type) — i.e. B_f is a
**building block with extra endomorphisms**, hence a point on the H₁₂ **badlocus**
(the "extra-endomorphism / rank-increase" component
`f⁶−f⁴−18ef²+27e²+16e`). Its reductions land on the badlocus at *every* prime, so
they are excluded from every M_ℓ; the generic curves in M_ℓ merely *share* B_f's
Frobenius trace (up to twist) without being B_f. This is the structural obstruction,
consistent with all the failed reconstructions.

**Status (final for this route).** The generic H₁₂ model provably cannot produce
B_f. Constructing it requires parametrizing the **badlocus itself** — the
1-dimensional locus of RM-by-√3 building blocks (a QM/Shimura-curve-type moduli
problem), or building B_f directly as a Q-surface from the elliptic building block
its inner twist descends to. That is a genuinely different (and lower-dimensional)
construction. NOTE: this building-block structure is itself the reason the earlier
"field of moduli Q / rational (e,f)" heuristic held.

## Artifacts
- `goal2_mod5_idx5.m` / `goal2_mod5_idx5_out.txt` — Step 1 + GRH pass (ℓ<5003 before timeout).
- `goal2_mod5_idx5_cert2.m` / `goal2_mod5_idx5_cert2_out.txt` — GRH pass [5000,12000].
- `goal2_mod5_Bf.m` / `goal2_mod5_Bf_out.txt` — B_f data: K_f=Q(√3), LMFDB `2.2.8.1-5000.1-j`, L-poly fingerprint.
- `bfs_family.m` — the validated BFS RM-by-√3 family (C_{1,2,−1} = sample check).
- `goal2_mod5_search.m` — blind box moduli search over Z[√2] (|trace| filter).
- `goal2_mod5_search2.m` — {2,5}-S-unit moduli search (c=1).
- `search_validate.m` — validates the fast point-counter against Magma.
- `goal2_mod5_igusaset.m` — per-prime P²(F_ℓ) scan: the decisive "not a quad twist of the family" diagnostic (0 matches at ℓ=7,23,31; 32 at ℓ=17).
- `goal2_mod5_h12.m` — Elkies–Kumar H₁₂ per-prime confirmation (matches at ALL of 7,17,23,31).
- `goal2_mod5_targets.m` / `goal2_mod5_targets_out.txt` — B_f (s,n) fingerprint at 20 split primes.
- `goal2_mod5_mcount.m` — |M_ℓ| counts (find small-|M| primes: 7,23 → 2).
- `goal2_mod5_crtQ.m` — rational CRT reconstruction of the H₁₂ moduli point (∏|M| method, height reach ~10⁵; blocked, see above).
- `goal2_mod5_curve.m` — build/verify a curve from a candidate (e,f).
- `goal2_mod5_loo2.m` — FAST builder (trace-prune) + leave-k-out + fresh-prime scoring; the run proving B_f is a building block on the H₁₂ badlocus (best score 3/7).
