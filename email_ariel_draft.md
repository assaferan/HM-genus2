Subject: HMF matching — swept the list (25 matched) + a question on the lifting step

Hi Ariel,

I went through your `sorted_output.jsonl` file. Summary below; full details and a
results table are in `results_goal1.md`.

## Candidate HMFs — 25 of 75 matched (all validated)

For each entry I match `A[p]` to `Ind_{G_F}^{G_Q}(ρ̄_f)` by finding a weight-[2,2]
Hilbert eigenform over F whose reduction reproduces the induced Frobenius data.
Two things were needed to make it work:

- **Level.** Your `serre_conductor_2d` is the prime-to-p conductor. For every
  mod-5 entry, 5 divides the conductor and is inert in Q(√2)/Q(√3), so the true
  newform level is Serre × 5⁴ (the 5-part comes back with exponent 2 and can't be
  traded for weight). The weight is forced to [2,2] by det = χ_p. For the mod-7
  entries 7 is a good prime, so the level is exactly the Serre norm.
- **Method.** Newform decomposition OOMs quickly, so I match instead by
  intersecting `ker(T_𝔭 − c_ℓ) mod p` over inert primes ℓ, with
  `c_ℓ = −[T² coeff of L_A at ℓ]` (skipping decomposition entirely). Each match
  is checked against a wrong-fingerprint control that must collapse to 0.

Matched: all mod-5 with real level norm ≤ ~405000 (idx 0–19, 22, 23) and the
three small mod-7 (idx 33, 53, 54). The other 50 have Hilbert spaces of dimension
62k–2.25M — beyond dense mod-p linear algebra. No shortcut removes the ×5⁴, so
pushing further needs a custom sparse mod-p Hecke implementation (Dembélé–Voight
directly over F_p). Happy to build that if the extra entries are worth it.

## Your mod-7 example (idx 33), in detail

Confirmed strongly: the full degree-4 characteristic polynomial of A mod 7 agrees
with `Ind ρ̄_f` at every good prime I tested (split and inert). One correction to
our plan, though: the matching form is **not rational** — its Hecke field has
degree 18. So the associated `B_f` is 18-dimensional, not an elliptic curve, and
the "construct `E_f`, form `B = Res(E_f)`, compare A[7] with B[7]" route doesn't
apply here. The eigenform gives the induced Frobenius data directly, which is all
we need for the comparison. (The small mod-5 forms, by contrast, have degree-2/4
Hecke fields — those `B_f` are genuinely low-dimensional if we want an explicit
construction somewhere.)

## Goal 2 — a question on the lifting step

I started on the modularity proof for idx 33 and hit exactly the subtlety you
flagged. From A's mod-7 Frobenius data I can show:
- the induced structure is airtight (0/224 inert primes fail the trace-0 test);
- `A[7]` is irreducible (`σ ≇ σ^τ` — the split-prime σ-traces are conjugate pairs
  in F_49 \ F_7);
- **σ has big image**: by Dickson it's not dihedral, not exceptional, and contains
  a transvection, so the projective image contains PSL₂(F₇).

So σ is a generic, big-image modular residual representation — good. **But** the
image of `A[7]` in GSp₄ is imprimitive (induced), so the standard BCGP
automorphy-lifting hypotheses (residual image containing Sp₄(F₇)) don't apply
off the shelf. My reading is that the right route for the induced case is to
**descend to GL₂/F** — lift σ via Hilbert modularity lifting (big image is fine
there), then induce the automorphic form back to GSp₄/Q by functoriality — rather
than lifting on the GSp₄ side. Does that match what you had in mind, or is there a
specific GSp₄ lifting statement you intended that allows imprimitive residual
image? The answer decides whether the Faltings–Serre step I'm setting up should
target the 4-dim GSp₄ isomorphism or the cleaner 2-dim `σ_A ≅ ρ̄_f` over F.

I'm building the 2-dim Faltings–Serre certificate over F in the meantime.

Best,
Eran
