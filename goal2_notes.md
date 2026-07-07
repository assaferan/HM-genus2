# Goal 2 — proving modularity of A (idx 33), roadmap + progress

**Target theorem.** `A = Jac(C)`, `C: y² = −10x⁵+15x⁴−30x³+10x²+60x+9`, is
modular (its L-function is automorphic / it matches a paramodular–GSp₄ object).

**Strategy.** Prove the residual isomorphism `A[7] ≅ Ind_{G_F}^{G_Q}(ρ̄_{f,λ})`
(`F=Q(√10)`, `f` = orbit-10 Hilbert newform, level norm 1500, weight [2,2], λ|7,
residue field F₄₉), then invoke automorphy lifting.

Ramification set `S = {2,3,5,7}` (`cond A = 2⁸·3·5⁵`, `+7`).

---

## Step 1 — the mod-7 image (DONE, `goal2_image.m`)

From A's char polynomial of Frobenius mod 7 (`goal2_image_out.txt`):

1. **Induced structure confirmed.** Every one of 224 primes inert in F (up to
   3000) has `a₁(A) ≡ 0 (mod 7)` — the anti-diagonal/swap-coset signature of an
   induced representation. 0 failures. At split primes `P_A = (T²−t₁T+ℓ)(T²−t₂T+ℓ)`.

2. **A[7] is irreducible.** At split ℓ the two σ-traces `t₁,t₂` are roots of an
   F₇-quadratic and land in **F₄₉∖F₇ as a conjugate pair** (`t₂=t₁⁷≠t₁`), so
   `σ ≇ σ^τ`; hence `Ind σ` is irreducible (4-dim, absolutely irreducible over F₇).

3. **σ has big image (Dickson).** For `σ: G_F → GL₂(F₄₉)`:
   - not **dihedral**: `tr = 0` at only 8/404 primes (≈2%; dihedral ⇒ ≈50%);
   - not **exceptional** (A₄/S₄/A₅): 24 distinct projective traces `t²/ℓ`;
   - a **transvection is present** (`t²/ℓ = 4`, an order-7 unipotent);
   ⇒ projective image **⊇ PSL₂(F₇)**. σ is absolutely irreducible with big image.

**Consequence.** σ is a "generic" modular residual representation — the worst
small-image sub-cases are excluded. This is exactly the favorable input for
modularity-lifting *over F* (Hilbert modular).

---

## Step 2 — the residual isomorphism as a THEOREM (open)

Currently `A[7] ≅ Ind ρ̄_f` is established as very strong numerical evidence:
42/42 split-prime traces + full degree-4 char polynomials mod 7 at ~30 primes
(`weilres_qsqrt10_c.m`), and the induced/image structure above.

To upgrade to a proof: **Faltings–Serre–Livné** for the pair of 4-dim mod-7
representations `ρ₁ = A[7]`, `ρ₂ = Ind ρ̄_f`, both unramified outside S with equal
determinant `χ₇`. Because the common image is a known finite group (imprimitive,
sitting in the normalizer of `GL₂(F₄₉)` inside `GSp₄(F₇)`), agreement of Frobenius
characteristic polynomials on an **effective** finite set of primes T — one whose
Frobenii meet every conjugacy class of the image — forces `ρ₁ ≅ ρ₂`. Concretely:
  (i) determine the image group G ⊆ GSp₄(F₇) and its conjugacy classes;
  (ii) exhibit primes in T covering every class, with matching char polys.
This is a finite, in-principle-mechanical certification; the remaining work is to
enumerate G's classes and pick the covering primes rigorously.

*Cleaner equivalent (the chosen route):* reduce to a **2-dimensional** comparison
over F. Since both sides restrict over `G_F` to `σ ⊕ σ^τ` resp. `ρ̄_f ⊕ ρ̄_f^τ`, it
suffices to prove `σ_A ≅ ρ̄_f` as 2-dim mod-7 reps of `G_F`. Lower-dimensional and
matches the big-image data from Step 1.

**Step 2 progress (DONE empirically, `goal2_step2.m`).** `σ_A` and `ρ̄_f` have
identical trace functions at every prime of F tested:
- **34/34 split primes** — `{t₁,t₂}` from A (roots of `Z²−a₁Z+(a₂−2ℓ)`) equals
  `{a_𝔭(f), a_𝔭'(f)} mod λ` as multisets in F₄₉;
- **5/5 inert primes** (ℓ = 11,17,19,23,29) — the unambiguous single-valued test
  `tr σ_A(Frob_{(ℓ)}) = −a₂(A,ℓ) = a_{(ℓ)}(f) mod λ` holds exactly in F₄₉.
Both sides have det = χ₇, and σ_A has big image (Step 1). This is the full
**input** to a Faltings–Serre certificate.

## Step 3 — precise images + the Goursat/Faltings–Serre reduction (DONE, `goal2_step3.m`)

**Both images are ⊇ SL₂(F₄₉) (rigorous).** From 155 Frobenius char polys of σ_A
(and 100 of ρ̄_f):
- **irreducible**: σ_A has 50 *nonsplit-semisimple* Frobenii (char poly irreducible
  over F₄₉), so its image lies in no Borel — σ_A is irreducible;
- **transvection**: an order-7 unipotent element is present;
- **trace field = F₄₉** (some traces lie in F₄₉∖F₇).
By **Dickson's classification**, irreducible + transvection + trace field F₄₉ ⟹
`im(σ_A) ⊇ SL₂(F₄₉)`. Same trace field for ρ̄_f; det = χ₇ surjects onto F₇*. So
both reps surject onto the same big group `G ⊇ SL₂(F₄₉)`.

**The reduction (Goursat + big image).** Let `H = im(σ_A × ρ̄_f) ⊆ G × G`; both
projections are onto G. Since `SL₂(F₄₉)` is quasi-simple, Goursat's lemma gives a
dichotomy:
1. **H is the graph of an automorphism** φ: G → G. Then `tr g = tr φ(g)` for all g.
   Aut(G) = inner · diagonal (PGL₂, trace-preserving) · ⟨Frobenius x↦x⁷⟩. The
   Frobenius part would force `tr g = (tr g)⁷` for all g, impossible since traces
   fill F₄₉ ⊋ F₇ — and it is killed concretely by trace agreement at any prime
   with trace in F₄₉∖F₇ (Step 2 has these). So φ is inner ⟹ **σ_A ≅ ρ̄_f**.
2. **H ⊇ SL₂(F₄₉) × SL₂(F₄₉)** ("independent"). Then H contains pairs (g,h) with
   `tr g ≠ tr h`, in fact `tr σ_A ≠ tr ρ̄_f` on a set of primes of density
   ≈ 1 − 1/49. Our verified agreement (Step 2: 39 primes both-sided, 0
   disagreements; would-be disagreement probability ≈ (1/49)³⁹ ≈ 10⁻⁶⁶ under
   equidistribution) excludes this case.

**Status of rigor.** The group theory (1)/(2) is a clean, correct reduction:
`σ_A ≅ ρ̄_f` follows once case (2) is excluded, i.e. once σ_A, ρ̄_f agree at *all*
primes. Turning "0 disagreements in 39 primes" into a theorem is the
effective-Chebotarev step, and here big image cuts both ways:

- **Livné's enumeration does NOT apply.** Livné is an ℓ=2 method: its deviation
  group is the elementary-abelian `Gal(F(√{S-units})/F)`. Our comparison is mod 7
  with big non-abelian image, and (both reps being residual) Faltings–Serre
  reduces to covering the conjugacy classes of `Gal(L/F)`, `L` = field cut out by
  `σ_A × ρ̄_f`. With images ⊇ SL₂(F₄₉), `|Gal(L/F)| ≤ |G|² ≈ 10¹¹`, so `L` is not
  explicitly constructible, and "generating the group from Frobenii" is
  insufficient (the agreement locus is not a subgroup).
- **GRH-conditional** effective Chebotarev gives a practical bound
  (`~(log disc)²`); **unconditional** effective Chebotarev (Lagarias–Odlyzko) is
  polynomial in `disc(L)` — astronomically large here. So a practical unconditional
  certificate via raw Galois comparison is out of reach for this big-image case.

**Upshot.** The residual iso is *morally* proved (clean reduction + `(1/49)³⁹`
odds) and can be made a **GRH-conditional theorem**. For genuine *unconditional*
modularity, don't force the raw comparison: use **automorphy lifting over F**
(Step 3 of strategy) — σ modular over F (big image is an asset for R=T; no
effective Chebotarev needed) then induce to GSp₄/Q. Big image obstructs the
unconditional Galois-theoretic route but *helps* the modularity route.

**Then:** automorphy lifting for σ over F (Hilbert modular, big image OK) and
induction back to GSp₄/Q — pending Ariel's confirmation of the route (see Step 3
of the strategy / the email).

---

## Step 3 — automorphy lifting (the real open question)

**The hurdle.** Although σ has big image in GL₂, the image of `A[7]` in `GSp₄(F₇)`
is **imprimitive** (induced by construction). The standard GSp₄ automorphy-lifting
theorems (BCGP and successors) assume the residual image is "enormous" / contains
`Sp₄(F_p)` — which an induced representation never satisfies. So BCGP does **not**
apply off-the-shelf here.

Two possible resolutions (for Ariel):
- **Descend to GL₂/F.** `A[7]=Ind σ` with σ modular (`= ρ̄_f`) and big-image; use
  the (very well-developed) **Hilbert modular lifting** for σ over F, then induce
  the resulting automorphic form back to GSp₄/Q by functoriality (automorphic
  induction / theta from GL₂/F to GSp₄/Q). This routes around the GSp₄ big-image
  hypothesis entirely and is the natural path for the induced case.
- Confirm whether a BCGP-type statement in the literature explicitly covers the
  residually-induced (imprimitive) case.

Note `A` is generic (`End=Z`), so its 7-adic image is open in `GSp₄(Z₇)`
(big); only the mod-7 *reduction* is accidentally imprimitive. So this is a
residual-image issue, not a geometric one.

**Question to Ariel:** which lifting statement did you intend for the induced
case — the GL₂/F-descend-then-induce route, or a specific GSp₄ theorem that
allows imprimitive residual image?

---

## Artifacts
- `goal2_image.m` / `goal2_image_out.txt` — Step 1 (image + induced structure).
- `weilres_qsqrt10_c.m` — the char-polynomial residual match (Step 2 evidence).
