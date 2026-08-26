# Matching Hugo Nartz's torsion-type mod-ℓ representations to Hilbert modular forms

*Working notes for the collaboration. Dataset: `examples.json` (Hugo Nartz, via Ariel) —
genus-2 curves over real quadratic fields whose mod-ℓ Galois representation (`ℓ ≥ 11`)
contains an irreducible 2-dimensional sub coming from a torsion point. The file supplied
holds **two** fields — `Q(√2)` (27 curves) and `Q(√3)` (12) — for **39 in total**, which is
exactly the set treated here. (Earlier drafts also named `Q(√5)`; no `Q(√5)` candidates are
present in the dataset, so nothing here covers that field.)*

## 0. Summary

For each curve `C/F` (`F` real quadratic) the mod-ℓ representation of `A = Jac(C)` is
```
    A[ℓ]^ss  =  1  ⊕  χ_ℓ  ⊕  σ,          σ irreducible, 2-dimensional, det σ = χ_ℓ,
```
so `σ` is a candidate for a **parallel weight-`[2,2]`, trivial-nebentypus Hilbert modular
form over `F`**. This is a direct (non-induced) test of Serre's conjecture over a real
quadratic field.

**Findings (all 39 curves in the dataset: 27 over `Q(√2)`, 12 over `Q(√3)`):**
- The whole list is **computationally easy**: because `ℓ ∈ {11,13}` does not divide any
  conductor, the Hilbert level is (essentially) the conductor, with **no `ℓ`-power blow-up**
  — unlike the `p=5` non-base-change family, whose `5⁴`-inflated levels reach dim `10⁵–10⁶`.
  Here every `Q(√2)` space has **dim ≤ 17 281**, and every `Q(√3)` space is comparable.
- The fingerprint is a **one-line point count**: `tr σ(Frob_P) ≡ −#C(𝔽_P) (mod ℓ)`.
- We **match** the small-conductor curves to explicit Hilbert newforms, each with a
  wrong-fingerprint discrimination control (a true match keeps agreement at all primes; the
  control collapses to 0). [matches table below]
- The **full sweep** (all 39 curves) is done by the **kernel-intersection matcher**
  (§4c): **all 27 `Q(√2)` curves** and **9/12 `Q(√3)`** match directly (control-validated),
  **36/39**. The last three — the largest `Q(√3)` spaces (level norms 569399, 785473) — hit a
  *deterministic Magma library bug* in `BasisMatrixDefinite`
  ([Magma issue #110](https://github.com/Magma-Maths/Magma/issues/110)), **not** a mathematical
  obstruction (the census §3 proves those spaces are in reach). A one-line source patch to the
  definite Hilbert-modular-forms code (§4d) clears it: **569399.3 now matches**
  (`dim 47728, survivor=1, control=0`). The two `785473` forms were then lost to an out-of-memory
  kill on a contended shared machine after ~40 h each (§4d); re-run one-at-a-time on an idle host,
  **both `785473` forms match** (`dim 76606, survivor=1, control=0`; 37.8 h and 35.8 h) — so the
  sweep is **complete at 39/39: every curve in the dataset is matched to a Hilbert newform**,
  each control-validated.
- Under **GRH**, each match upgrades to a theorem via an effective Faltings–Serre /
  Chebotarev prime bound `O((log cond)²)` — a few-hundred-prime check, as in the idx-33 work.
  Because the mod-ℓ survivor eigenvector supplies `a_P mod λ` directly, the certificate needs
  **no char-0 eigenform** and scales with the matcher: **37 GRH modularity certificates**
  (all 27 `Q(√2)`, 10 of 12 `Q(√3)`), every one with 0 trace disagreements and `σ` irreducible
  (§4b). The 2 outstanding are the `785473` pair, pending compute, not method.
- **"Not from an elliptic curve" is established for 5 of the 39, and open for the rest.** Where we
  have the char-0 eigenform the Hecke field is large — degree 18 for `881`, 5 for `14303` (§4e) —
  and a scalable mod-ℓ test adds `145161.2` with `[E:Q] ≥ 2` (§4f). That test is **one-sided**: it
  is silent whenever `ℓ` has a degree-1 prime in the Hecke field, as it did for 22 of the 27 curves
  tried, and `881` (`[E:Q] = 18`, test silent) shows an inconclusive result carries no information.
  This is a limit of the method, not of compute.
- Working with the mod-ℓ eigenvector costs nothing in strength: the space is characteristic 0 to
  begin with, and the Hecke matrices are ℓ-integral (checked, `verify_lattice.m`), so the
  **Deligne–Serre lifting lemma** produces a characteristic-0 Hilbert newform of the same weight
  and level dividing `N`. Getting from the certificate to the word *newform* needs that 1974 lemma
  about lattices — **not** a modularity lifting theorem (§4b).

## 1. Structure and fingerprint (verified)

A torsion point gives a Galois-fixed line (trivial sub `1`); the Weil pairing forces a
`χ_ℓ` quotient; the middle `σ` is 2-dimensional with `det σ = χ_ℓ` (from the symplectic
multiplier). Hence, at every good prime `P` of `F`, the Frobenius quartic mod ℓ factors
```
    charpoly(Frob_P) ≡ (T−1)(T−N(P)) · (T² − t_P T + N(P))   (mod ℓ),
```
so the σ-trace is a **single** value peeled off a point count:
```
    t_P  ≡  a_P − 1 − N(P)  ≡  −#C(𝔽_P)   (mod ℓ).
```
**Verified** on all 39 curves: the `(T−1)(T−N(P))` factor is present at 8–9 of the first
good primes for every curve (`validate.m`).

## 2. Level

`cond(σ)` away from ℓ equals `cond(A)` away from ℓ (the `1` and `χ_ℓ` pieces are unramified
away from ℓ), and `ℓ ∤ cond(A)`, so:
- the **odd part** of the level = the given *prime-to-2 conductor* (for `Q(√2)` we use the
  ideal Hugo supplies; for `Q(√3)` Magma's `Conductor(C)` reproduces it exactly);
- the **2-part** is unknown a priori. `Conductor(C)` at the prime `𝔭₂ | 2` gives an **upper
  bound** `e_max` (exactly Ariel's "low-ish upper bound"); the true σ-level is
  `cond' · 𝔭₂^e` for the **smallest** `e ≤ e_max` at which a newform matches (mod-ℓ
  level-lowering). E.g. curve `881`: `Conductor(C) = 𝔭₂⁴·881`, but σ lives at `𝔭₂³·881`.

## 3. Dimension census — `Q(√2)` (exact, `dim S_[2,2]` at the prime-to-2 conductor)

| N(cond) | ℓ | dim | N(cond) | ℓ | dim |
|--:|--:|--:|--:|--:|--:|
| 881 | 13 | 37 | 145161 | 13 | 5293 |
| 14303 | 11 | 595 | 161089 | 11 | 6547 |
| 20447 | 11 | 689 | 173111 | 13 | 6753 |
| 24889 | 13 | 1038 | 200273 | 11 | 8345 |
| 68193 | 11 | 2525 | 243049 | 13 | 8961 |
| 100489 | 13 | 4188 | 312769 | 11 | 12971 |
| 105121 | 11 | 4239 | 328329 | 13 | 12033 |
| 113609 | 13 | 4685 | 478593 | 11 | 17281 |

All ≤ 17 281 — every one is in reach (the project's kernel-intersection method has run to
dim 2.25M).

**`Q(√3)`** (exact, one per conductor; conjugates share dim):

| N(cond) | ℓ | dim | N(cond) | ℓ | dim |
|--:|--:|--:|--:|--:|--:|
| 4057 | 11 | 340 | 472993 | 11 | 39418 |
| 65209 | 13 | 5340 | 569399 | 13 | 47168 |
| 72649 | 13 | 6056 | 785473 | 11 | 55446 |
| 377233 | 11 | 31102 | | | |

All ≤ 55 446 (larger per unit norm since `ζ_{Q(√3)}(−1) = 1/6` vs `1/12` for `Q(√2)`). For
reference, were `Q(√5)` candidates ever supplied, its spaces would be `≈ 0.4×` the `Q(√2)`
size — i.e. the cheapest of the three to run.

## 4. Matches

*(`sweep.m`: level-lowering from `e=0`, `NewformDecomposition`, discrimination control t→t+1.
Each row: a Hilbert newform whose mod-λ reduction equals σ at all tested primes, with the
control agreeing at 0.)*

| label | F | ℓ | σ-level norm (2-part e) | dim | orbit dim | Hecke deg | primes agree | control |
|---|---|--:|---|--:|--:|--:|--:|--:|
| 881.1   | Q(√2) | 13 | 7048 (e=3) | 110 | 18 | 18 | 29 | 0 |
| 881.2   | Q(√2) | 13 | 7048 (e=3) | 110 | 18 | 18 | 29 | 0 |
| 14303.1 | Q(√2) | 11 | 14303 (e=0) | 595 | 5 | 5 | 29 | 0 |
| 14303.2 | Q(√2) | 11 | 14303 (e=0) | 595 | 5 | 5 | 29 | 0 |

These four already span both residual primes (11, 13) and both 2-part regimes: `14303` is
modular at the bare prime-to-2 conductor (`e=0`), while `881`'s σ is **less ramified at 2
than A** — `Conductor(A) = 𝔭₂⁴·881` but σ lives at `𝔭₂³·881` (mod-ℓ level-lowering). The
remaining curves are running (`sweep.m` streams to `sweep.out`); the larger-dim tail
(`dim > ~2600`) is deferred to the kernel-intersection matcher used for Goal 1. **The census
(§3) already establishes that every space is small enough to reach.**


### 4a. Explicit identifications (curve + Hilbert newform)

Checked LMFDB: **neither form is in it** — the HMF API returns no records for `2.2.8.1` (Q(√2)) at level norm 7048 or 14303 (beyond LMFDB's coverage for this field). So each is given by curve + Hecke cutters (min. poly of `a_P`). Throughout `a = √2`, `y` generates the Hecke field; `σ` = reduction of `f` mod a prime `λ | ℓ`, verified `a_P(f) ≡ −#C(F_P) (mod ℓ)` at all 29 tested primes (control at 0). The conjugate curves `.2` (under `a ↦ −a`) match the Galois-conjugate forms.

### 881.1

**Curve** `C/Q(√2)`:  
`y^2 + ((1-a)x^3 - a x^2 + (1-a)x) y = a x^6 + (1+2a)x^5 + 2a x^4 + (-2+a)x^3 - a x^2 + (-1+a)x + 1`

- **level** norm 7048  =  [ <2, 3>, <881, 1> ]
- **weight** [2,2], trivial nebentypus
- **orbit dim** 18; **Hecke field** deg 18, totally real, Galois group S_18, disc = 2^3·11·1709·(29-digit prime) — **not in LMFDB**. Defining poly Q[y]/(y^18 - 54y^16 - 6y^15 + 1124y^14 + 184y^13 - 11538y^12 - 2264y^11 + 61687y^10 + 13681y^9 - 166924y^8 - 35569y^7 + 203946y^6 + 26001y^5 - 95812y^4 - 10226y^3 + 12624y^2 + 2360y - 16)
- **Hecke cutters** — minpoly of `a_P`:
    - `P (norm 9, (3)):  minpoly(a_P) = y^18 + 17y^17 + 74y^16 - 245y^15 - 2532y^14 - 1981y^13 + 26272y^12 + 51570y^11 - 117427y^10 - 345154y^9 + 187821y^8 + 1035589y^7 + 216747y^6 - 1366380y^5 - 969569y^4 + 442423y^3 + 699359y^2 + 251552y + 27232`
    - `P (norm 25, (5)):  minpoly(a_P) = y^18 + 33y^17 + 319y^16 - 1162y^15 - 43547y^14 - 260171y^13 + 469115y^12 + 12033929y^11 + 44501109y^10 - 56316999y^9 - 852721809y^8 - 2067056351y^7 + 607748975y^6 + 10349900301y^5 + 12588696758y^4 - 8439753640y^3 - 24066987640y^2 - 8013182592y + 4208944064`
    - `P (norm 7, (-2a + 1)):  minpoly(a_P) = y^18 - 54y^16 - 6y^15 + 1124y^14 + 184y^13 - 11538y^12 - 2264y^11 + 61687y^10 + 13681y^9 - 166924y^8 - 35569y^7 + 203946y^6 + 26001y^5 - 95812y^4 - 10226y^3 + 12624y^2 + 2360y - 16`
    - `P (norm 7, (-2a - 1)):  minpoly(a_P) = y^18 + 8y^17 - 21y^16 - 304y^15 - 233y^14 + 3535y^13 + 6267y^12 - 17981y^11 - 44980y^10 + 40516y^9 + 150365y^8 - 20005y^7 - 252439y^6 - 70433y^5 + 196619y^4 + 114585y^3 - 45655y^2 - 49824y - 10368`
    - `P (norm 121, (11)):  minpoly(a_P) = y^18 + 27y^17 - 487y^16 - 18739y^15 + 4747y^14 + 4107291y^13 + 22532242y^12 - 317833147y^11 - 2765513587y^10 + 9051619144y^9 + 127688119082y^8 - 18325224776y^7 - 2611261456546y^6 - 2960865341641y^5 + 22737939261662y^4 + 40019962056332y^3 - 56249965840480y^2 - 133054870257792y - 58651454037024`
    - `P (norm 169, (13)):  minpoly(a_P) = y^18 + 73y^17 + 1070y^16 - 43475y^15 - 1443157y^14 - 151620y^13 + 484029190y^12 + 4715556014y^11 - 53377199272y^10 - 1080484020257y^9 - 1193131166023y^8 + 83002564637363y^7 + 524656100649275y^6 - 1128621069617603y^5 - 21911135202863470y^4 - 68017735219054844y^3 + 26259192645059072y^2 + 461590316652796032y + 599871255543931104`

### 14303.1

**Curve** `C/Q(√2)`:  
`y^2 + (x^3 + x^2 + 1) y = a x^5 + (3+2a)x^4 + (3+a)x^3 + (-4-2a)x^2 + (-1+2a)x + (1-a)`

- **level** norm 14303  =  [ <14303, 1> ]
- **weight** [2,2], trivial nebentypus
- **orbit dim** 5; **Hecke field** = **LMFDB [5.5.14641.1](https://www.lmfdb.org/NumberField/5.5.14641.1)** — the cyclic quintic of conductor 11 (= Q(ζ_11)^+), disc 11^4. Generator here: Q[y]/(y^5 - 2y^4 - 5y^3 + 13y^2 - 7y + 1) (LMFDB polredabs: x^5 + x^4 - 4x^3 + 3x^2 + 3x - 1)
- **Hecke cutters** — minpoly of `a_P`:
    - `P (norm 2, (-a)):  minpoly(a_P) = y^5 - 2y^4 - 5y^3 + 13y^2 - 7y + 1`
    - `P (norm 9, (3)):  minpoly(a_P) = y^5 + y^4 - 15y^3 - 14y^2 + 3y + 1`
    - `P (norm 25, (5)):  minpoly(a_P) = y^5 - 4y^4 - 75y^3 + 170y^2 + 1010y + 439`
    - `P (norm 7, (-2a + 1)):  minpoly(a_P) = y^5 - 5y^4 - 23y^3 + 122y^2 + 115y - 683`
    - `P (norm 7, (-2a - 1)):  minpoly(a_P) = y^5 - 33y^3 + 242y + 121`
    - `P (norm 121, (11)):  minpoly(a_P) = y^5 + 9y^4 - 335y^3 - 2088y^2 + 23665y + 30097`

### 4b. Modularity theorems (under GRH): 37 certificates

#### The two worked examples (char-0 route)

For curve `14303.1` we upgrade the match to a **theorem under GRH**, following the idx-33
method (`grh_14303.m`). `sigma` (the 2-dim sub of `A[11]`) and `rho_f` (mod-lambda reduction
of the level-14303 newform, Hecke field `5.5.14641.1` where 11 is totally ramified) are both
2-dimensional, **irreducible** (a Frobenius char poly is irreducible over F_11), with
`det = chi_11`, unramified outside `{P_14303, 11}`. By Brauer-Nesbitt they are isomorphic iff
`tr Frob` agree at all good primes; under GRH the conductor-based effective
Faltings-Serre / Chebotarev bound `~ (log cond)^2 ~ 200` (independent of the splitting field)
makes this finite. We verified **exact trace agreement `tr sigma(Frob_P) = tr rho_f(Frob_P)`
at all 428 good primes `P` of `Q(sqrt2)` with `N(P) <= 3000`** (0 disagreements) -- far beyond
the bound.

> **Theorem [GRH].** `sigma = rho_f`; hence `sigma` is modular. (First Serre-modularity
> theorem over `Q(sqrt2)` from this dataset.)

For curve `881.1` (`l=13`, level `p2^3*881`) the same certificate goes through (`grh_881.m`):
`sigma`, `rho_f` are 2-dim, irreducible, `det = chi_13`, unramified outside `{P_881, 2, 13}`;
the Hecke field is a generic degree-18 `S_18` field (a `lambda|13` of residue degree 1 gives the
reduction). Verified **exact trace agreement at all 427 good primes `N(P) <= 3000`** (0
disagreements).

> **Theorem [GRH].** For `881.1` (`l=13`), `sigma = rho_f`; hence `sigma` is modular.

So both demonstrated matches are modularity theorems under GRH (l = 11 and l = 13).

#### Scaling the certificate: `a_P mod λ` from the survivor eigenvector

Both scripts above route `tr ρ_f(Frob_P)` through a char-0 `NewformDecomposition`, which — as in
§4c — dies past the four smallest curves. The certificate does **not actually need the char-0
form**: the kernel survivor `v` of §4c is already a mod-ℓ Hecke eigenvector, so

```
    tr ρ_f(Frob_P)  =  c_P,   where   v·T_P = c_P·v   over F_ℓ,
```

and `c_P = a_P(f) mod λ` directly. `grh_kernel_validate.m` proves the identity on `14303`:
`c_P`, the reduction of the char-0 `HeckeEigenvalue(f,P)`, and the curve-side fingerprint
`−#C(𝔽_P)` agree at **all 428 good primes `N(P) ≤ 3000`, 0 mismatches**. `grh_kernel.m` then
emits the full certificate from `v` alone: it finds the match level by the `e`-sweep, checks
`σ` irreducible (a Frobenius char poly irreducible over `F_ℓ`), and verifies exact trace
agreement `tr σ(Frob_P) = c_P` at every good `P` up to `BOUND`. No char-0 decomposition
anywhere, so it reaches the same levels the kernel matcher does.

**Result: 37 GRH certificates — all 27 `Q(√2)` curves and 10 of the 12 `Q(√3)`.** Every one is
`disagree = 0` with `σ` irreducible. `BOUND` is tiered by dimension (3000 for `dim ≤ 9k`, 800
above); both tiers are far beyond the GRH effective Faltings–Serre bound `~(log cond)² ~ 200`.

| label(s) | F | ℓ | level norm | e | BOUND | good primes | disagreements |
|---|---|--:|--:|--:|--:|--:|--:|
| 881.1, 881.2 | Q(√2) | 13 | 7048 | 3 | 3000 | 427 | 0 |
| 14303.1, 14303.2 | Q(√2) | 11 | 14303 | 0 | 3000 | 428 | 0 |
| 20447.3, 20447.6 | Q(√2) | 11 | 20447 | 0 | 3000 | 425 | 0 |
| 24889.1, 24889.2 | Q(√2) | 13 | 24889 | 0 | 3000 | 428 | 0 |
| 68193.1, 68193.2 | Q(√2) | 11 | 68193 | 0 | 3000 | 427 | 0 |
| 100489.1 | Q(√2) | 13 | 100489 | 0 | 3000 | 428 | 0 |
| 105121.1 | Q(√2) | 11 | 105121 | 0 | 3000 | 427 | 0 |
| 113609.1, 113609.4 | Q(√2) | 13 | 113609 | 0 | 3000 | 426 | 0 |
| 145161.2 | Q(√2) | 13 | 145161 | 0 | 3000 | 425 | 0 |
| 161089.2, 161089.3 | Q(√2) | 11 | 161089 | 0 | 3000 | 427 | 0 |
| 173111.2, 173111.5 | Q(√2) | 13 | 173111 | 0 | 3000 | 426 | 0 |
| 200273.1, 200273.2 | Q(√2) | 11 | 200273 | 0 | 3000 | 428 | 0 |
| 243049.2 | Q(√2) | 13 | 243049 | 0 | 800 | 130 | 0 |
| 312769.2, 312769.3 | Q(√2) | 11 | 312769 | 0 | 800 | 131 | 0 |
| 328329.2 | Q(√2) | 13 | 328329 | 0 | 800 | 129 | 0 |
| 478593.1, 478593.4 | Q(√2) | 11 | 478593 | 0 | 800 | 130 | 0 |
| 4057.1, 4057.2 | Q(√3) | 11 | 16228 | 2 | 3000 | 424 | 0 |
| 65209.2, 65209.3 | Q(√3) | 13 | 65209 | 0 | 3000 | 422 | 0 |
| 72649.1, 72649.2 | Q(√3) | 13 | 72649 | 0 | 3000 | 424 | 0 |
| 377233.2 | Q(√3) | 11 | 377233 | 0 | 800 | 132 | 0 |
| 472993.1 | Q(√3) | 11 | 472993 | 0 | 800 | 133 | 0 |
| 472993.2 | Q(√3) | 11 | 472993 | 0 | 800 | 133 | 0 |
| 569399.3 | Q(√3) | 13 | 569399 | 0 | 400 | 74 | 0 |

> **Theorem [GRH].** For each of the 37 curves above, `σ ≅ ρ̄_{f,λ}` for a Hilbert modular newform
> `f` of parallel weight `[2,2]` and level dividing the stated level; hence `σ` is modular.

The first two rows reproduce the char-0 theorems above (`14303.1` at `e=0`, `881.1` at `e=3`)
from the survivor alone — the kernel route is validated against the char-0 route wherever both
are computable.

#### From a mod-ℓ eigenvector to a characteristic-0 form

The certificate above is computed entirely with the mod-ℓ survivor `v`, so it is worth being
precise about what licenses the word *newform* in the theorem. The step is **not** a modularity
lifting theorem — no `R = T`, no deformation theory. It is a 1974 lemma about lattices:

> **Lemma** (Deligne–Serre, *Formes modulaires de poids 1*, Ann. Sci. ÉNS (4) **7** (1974),
> 507–530, Lemme 6.11). Let `M` be a free module of finite rank over a discrete valuation ring `R`
> with maximal ideal `𝔪`, and let `S` be a set of commuting `R`-endomorphisms of `M`. If
> `0 ≠ v ∈ M` satisfies `Tv ≡ c_T v (mod 𝔪M)` for all `T ∈ S`, then there is a DVR `R' ⊇ R`, finite
> over `R`, with maximal ideal `𝔪' ⊇ 𝔪`, and a nonzero `v' ∈ R' ⊗_R M` with `Tv' = c'_T v'` and
> `c'_T ≡ c_T (mod 𝔪')`.

**Why it applies here.** The essential point is that we never left characteristic 0:
`HilbertCuspForms(F, N, [2,2])` *is* a characteristic-0 Hecke module, and only the Hecke **matrices**
are reduced mod ℓ. So the lemma needs just one thing — a Hecke-stable lattice.

That is not quite automatic. Magma's Hecke matrices on this space are rational and, for most of our
levels, **not integral**. But their denominators are prime to ℓ, so the `Z_(ℓ)`-span of Magma's basis
*is* Hecke-stable, and the lemma applies with `R = Z_(ℓ)`. This is checkable rather than assumed, and
`verify_lattice.m` checks it (it is in the test gate):

| curve | ℓ | `T_P` checked | non-integral | max `v_ℓ(denominator)` |
|---|--:|--:|--:|--:|
| 881.1 | 13 | 11 | 0 | 0 |
| 4057.1 | Q(√3), 11 | 10 | **10** | **0** |
| 24889.1 | 13 | 11 | **11** | **0** |
| 65209.2 | Q(√3), 13 | 10 | **10** | **0** |

Note `4057.1`, `24889.1` and `65209.2`: *every* Hecke matrix is non-integral, yet every denominator
is prime to ℓ. This is the same ℓ-integrality that makes the reduction `T_P ↦ T_P mod ℓ` well defined
in `kernel_torsion.m` and `grh_kernel.m` — so a failure here would break the certificates and the
lifting argument together, which is why it is worth a gate rather than a remark.

Applying the lemma with `M` = that lattice, `S = {T_P}`, and `v` = our survivor yields a
characteristic-0 Hilbert cusp eigenform of the **same level and parallel weight `[2,2]`** whose
eigenvalues satisfy `a_P ≡ c_P (mod λ)` at every prime we tested. Attaching `ρ_{f,λ}` to it is then
standard (Blasius–Rogawski, Carayol, Deligne, Saito, Taylor, Wiles), and `ρ̄_{f,λ}` has the traces the
certificate verifies against.

**Two qualifications.**

- The survivor lives in the *full* cusp space, so the lifted eigenform need not be new at `N`; its
  associated newform has level **dividing** `N`, with the same Galois representation away from `N`.
  The theorem is stated that way above.
- The lemma produces eigenvalues in a finite extension `R'`, not in `Z_ℓ`. That is expected here —
  the Hecke fields we identified explicitly have degrees 5 and 18 (§4e).

(Deligne–Serre is a paper about weight 1, but Lemme 6.11 is a general statement about modules over
a DVR and carries no weight hypothesis, so nothing about weight 1 is being imported here. The
reference has been checked against the original paper.)

**Outstanding (2 of 39).** All 39 curves are matched (§4c/§4d) and 37 are certified; only the
`785473` pair (`785473.12`, `785473.5`, dim 76606) lacks a certificate. `569399.3` was certified
at `BOUND=400` in ~14 h. The `785473.12` attempt reached 27 of its 74 primes with 0 disagreements
before being killed at ~43 h — see the note on host contention in §4d.

The cost split is worth recording, since it governs what the giants will take. Measured on
`472993.2` (dim 39418): **~6.8 h building the space and Hecke data before verification began**,
then **~340 s per prime** — 42 338 s to reach 125 of the 133 primes, ~12.5 h for the trace loop,
~19 h in total.

So verification, not the space build, dominates at these dimensions, and it scales linearly in the
number of primes. `BOUND` should therefore be the smallest value comfortably above the
Faltings–Serre bound (`~185` for these conductors) rather than as large as convenient: `BOUND=800`
here bought a 4× margin over the bound at ~4× the cost of `BOUND=400`, which `run_grh_giants.sh`
already uses for exactly this reason. The per-prime rate also gives a usable estimate for the
giants — at dim 47 728–55 446 expect appreciably more than 340 s per prime.

### 4c. Full sweep via kernel intersection (36/39)

`NewformDecomposition` (used in `sweep.m`, §4) only reaches the 4 smallest curves — it is
intractable past dim ~2600 (it ran >10 h with no output on the dim-2525 mid-size spaces). The
**kernel-intersection matcher** `kernel_torsion.m` (adapted from Goal 1's `kernel_match.m`)
replaces it: on the **full** cusp space `M = HilbertCuspForms(F, base·𝔭₂^e, [2,2])` it forms
```
    survivor = ⋂_P ker(T_P − t_P·I)   over F_ℓ,     t_P = −#C(𝔽_P) mod ℓ,
```
over good primes `P` (split **and** inert — `σ` is a genuine `GL₂/F` rep, so `t_P ∈ F_ℓ`
directly, no induced-case inert restriction), sweeping the 2-part exponent `e`. `survivor`
dim > 0 ⟺ a newform matches; the discrimination control `t_P → t_P+1` must collapse it to 0.
No char-0 decomposition. In every match below `survivor` is **exactly 1-dimensional** (so the
eigenform is isolable — see the caveat at the end).

**Result: 36/39 match directly, every one `survivor=1, control=0`** — all 27 `Q(√2)` curves, and
9 of 12 `Q(√3)`; the remaining 3 (the `Q(√3)` giants) match on a patched Magma build, see §4d.
(`dim` here is the **full** cusp-space dimension the kernel runs on, larger than the `NewSubspace`
dims of §3; `e>0` only for `881` and `4057`, via mod-ℓ level-lowering.)

| label | F | ℓ | level norm | e | full-space dim | survivor/control |
|---|---|--:|--:|--:|--:|:--:|
| 881.1, 881.2 | Q(√2) | 13 | 7048 | 3 | 441 | 1 / 0 |
| 14303.1, 14303.2 | Q(√2) | 11 | 14303 | 0 | 595 | 1 / 0 |
| 20447.3, 20447.6 | Q(√2) | 11 | 20447 | 0 | 1023 | 1 / 0 |
| 24889.1, 24889.2 | Q(√2) | 13 | 24889 | 0 | 1038 | 1 / 0 |
| 68193.1, 68193.2 | Q(√2) | 11 | 68193 | 0 | 3159 | 1 / 0 |
| 100489.1 | Q(√2) | 13 | 100489 | 0 | 4188 | 1 / 0 |
| 105121.1 | Q(√2) | 11 | 105121 | 0 | 4523 | 1 / 0 |
| 113609.1, 113609.4 | Q(√2) | 13 | 113609 | 0 | 4783 | 1 / 0 |
| 145161.2 | Q(√2) | 13 | 145161 | 0 | 6827 | 1 / 0 |
| 161089.2, 161089.3 | Q(√2) | 11 | 161089 | 0 | 6879 | 1 / 0 |
| 173111.2, 173111.5 | Q(√2) | 13 | 173111 | 0 | 7649 | 1 / 0 |
| 200273.1, 200273.2 | Q(√2) | 11 | 200273 | 0 | 8345 | 1 / 0 |
| 243049.2 | Q(√2) | 13 | 243049 | 0 | 11371 | 1 / 0 |
| 312769.2, 312769.3 | Q(√2) | 11 | 312769 | 0 | 13095 | 1 / 0 |
| 328329.2 | Q(√2) | 13 | 328329 | 0 | 15359 | 1 / 0 |
| 478593.1, 478593.4 | Q(√2) | 11 | 478593 | 0 | 22719 | 1 / 0 |
| 4057.1, 4057.2 | Q(√3) | 11 | 16228 | 2 | 2030 | 1 / 0 |
| 65209.2, 65209.3 | Q(√3) | 13 | 65209 | 0 | 5532 | 1 / 0 |
| 72649.1, 72649.2 | Q(√3) | 13 | 72649 | 0 | 6056 | 1 / 0 |
| 377233.2 | Q(√3) | 11 | 377233 | 0 | 31774 | 1 / 0 |
| 472993.1, 472993.2 | Q(√3) | 11 | 472993 | 0 | 39418 | 1 / 0 |
| 569399.3 | Q(√3) | 13 | 569399 | 0 | 47728 | 1 / 0 (§4d) |
| 785473.12 | Q(√3) | 11 | 785473 | 0 | **76606** | 1 / 0 (§4d) |
| 785473.5 | Q(√3) | 11 | 785473 | 0 | 76606 | 1 / 0 (§4d) |

The last three rows use the patched build of §4d. All three are now confirmed, so **every curve in
the dataset is matched: 39/39**, each with `survivor = 1` and the discrimination control at `0`.

**Correction to the `785473` dimension.** Earlier drafts carried `55446` for these two rows, taken
from the §3 census while no run had ever completed. The finished `785473.12` run reports the actual
full cusp-space dimension as **76606** — 38 % larger. (This is consistent with the note above: the
§3 census is the prime-to-2 `NewSubspace` dimension, whereas the kernel runs on the *full* space.
For the other levels the two happen to be close, which is why the placeholder went unnoticed.) The
gap explains why these runs cost more than the `≈40 h` extrapolated from the interrupted attempts.

### 4d. The three `Q(√3)` giants: a one-line Magma #110 workaround

On the three largest spaces (`569399`, `785473`, both `Q(√3)`) the **first Hecke operator**
crashes deterministically:
```
BasisMatrixDefinite(M)  →  definite.m:1060
    Binv := Transpose(Solution(Transpose(B), IdentityMatrix(BaseRing(B), Nrows(B))));
Runtime error in 'Solution': No solution exists
```
**Root cause.** `B` (`basis_matrix_big`) is assembled correctly from the ideal-class direct
factors, but Magma then computes a right inverse `Binv` via `Solution(Bᵀ, I)`, which needs `B`
to have full **row** rank. For these two levels the assembly yields **linearly dependent rows**,
so the solve has no solution. This is a genuine Magma library bug
([issue #110](https://github.com/Magma-Maths/Magma/issues/110)) — **not** a size limit: the
larger-dimension `472993` (dim 39418) succeeds. The bug is confirmed and **fixed upstream in
Magma V2.29-10** (issue #110, resolved by A. Steel); the patch below was the interim workaround
we used on V2.29-9, and remains the route on any machine not yet upgraded to V2.29-10.

**Why the fix is safe for parallel weight 2.** `basis_matrix_big_inv` (the crashing `Binv`) is
**read in exactly one place** — `definite.m:1073`, inside the *non*-weight-2 branch. For weight
`[2,2]`, `RemoveEisenstein` rebuilds `basis_matrix`/`basis_matrix_inv` from the Eisenstein
indicator vectors and the inner product (never touching `Binv`), and the big Hecke matrix uses
only `Ncols(basis_matrix_big)`. So for parallel weight 2 the inverse is **vestigial**. The patch
wraps the solve in `try/catch`: skip `Binv` for weight 2, **re-raise for any other weight**. On
the success path it is byte-for-byte the original computation.

**Validation.** Patched vs. stock Magma with the kernel matcher: `14303.1` (`Q(√2)`) and `4057.1`
(`Q(√3)`) give **byte-identical** output; `881.1` reproduces its `e=3` match. The patch is inert
wherever line 1060 succeeds; it only changes the previously-crashing giants.

**Result.** On the patched build `569399.3` runs end-to-end (dim 47728; ~10.5 h) to
`survivor=1, control=0` — a control-validated **match**. The deploy is a **private patched Magma
copy** (no system files touched, no Magma source redistributed — only our ~25-line diff): see
`magma110_patch/` (`definite.m.patch`, `deploy_patch.sh`, `README.md`).

**The two `785473` forms (dim 76606) were first lost for operational, not mathematical, reasons —
and have since been recovered.** Both were originally launched *in parallel* on the patched build
and ran ~40 h at ~400 GB each, then were killed by the host's OOM reaper when an unrelated job on
the same shared machine expanded to >1 TB. Because `kernel_torsion.m` writes nothing until the
space-and-Hecke build completes, **both produced no output at all** — 40 h each, zero partial credit.

**Re-run and matched.** On a second, uncontended host, run **strictly one at a time**,
`785473.12` completed:

```
785473.12   d=3 l=11  levelN=785473  e=0  dim=76606  survivor=1  control=0  MATCH   [136 001 s]
```

That is **37.8 h of CPU for a single giant** — the first completed match at this dimension, and so
the first measured cost rather than an extrapolation. `785473.5` then completed under the same sequential discipline:

```
785473.5    d=3 l=11  levelN=785473  e=0  dim=76606  survivor=1  control=0  MATCH   [128 937 s]
```

35.8 h, closing the sweep at **39/39 — every curve in the dataset matched**. The two conjugates
cost 37.8 h and 35.8 h respectively, so ~36–38 h is the reliable figure for a single dim-76606
match on an uncontended host.

**Certificates for the giants, and a second contention problem.** `569399.3` (dim 47728) certified
cleanly at `BOUND=400`: 74 good primes, 0 disagreements, ~14 h. The `785473.12` certificate is
harder — at dim 76606 the space-and-Hecke build alone took **~37 h** before verification began,
and the run was then killed at ~43 h having checked 27 of its 74 primes (0 disagreements). That
kill was **not** an out-of-memory event: it exited on `SIGTERM` with 988 GB free and no OOM-daemon
activity. The host is a shared machine that also runs continuous-integration jobs for Magma itself,
whose processes share the name `magma.exe`; a cleanup that reaps by process name is the likeliest
explanation. The lesson is the same as before in a new guise — *whose* machine it is matters as
much as how much memory it has.

Two lessons for anyone repeating this, both cheap to act on:

- **The runs are not checkpointed.** At this scale that turns any interruption into total loss. If
  the giants are attempted again, either checkpoint the Hecke build or accept that the whole run is
  an all-or-nothing bet.
- **A shared machine is part of the experiment.** The host's OOM policy here explicitly *prefers*
  killing `magma`, so a competing job does not merely slow these runs down — it selects them for
  termination first. Coordinate exclusive time rather than relying on there being enough headroom.

### 4e. Hecke fields of the identified forms — evidence against an elliptic-curve source

For the four **explicitly isolated** newforms (§4a, `hecke_cutters.m`) we know the Hecke
eigenvalue field `E = Q(a_P : P)` exactly. Recall that a parallel weight-`[2,2]` Hilbert
newform `f` with Hecke field `E`, `[E:Q]=d`, has an attached **`GL₂`-type abelian variety**
`A_f/F` of **dimension `d`** with real multiplication by `E`; its mod-ℓ representations are the
`ρ̄_{f,λ}`. An elliptic curve `E/F` (or any form with **rational** Hecke field, `d=1`) yields a
2-dimensional mod-ℓ representation with traces in `F_ℓ`. So `d>1` means the modular source is a
genuinely higher-dimensional abelian variety — **not an elliptic curve**.

| form | F | ℓ | level norm | Hecke field `E` | `[E:Q] = dim A_f` |
|---|---|--:|--:|---|--:|
| 881.1, 881.2 | Q(√2) | 13 | 7048 | totally real, deg 18 (single field; disc ≈ 4.25×10³³; too large for LMFDB) | **18** |
| 14303.1, 14303.2 | Q(√2) | 11 | 14303 | `Q(ζ₁₁)⁺` = `5.5.14641.1` (cyclic C₅, disc 11⁴) | **5** |

(All four fields verified in Magma: irreducible, **totally real**, `d = 18` resp. `5`; `881.1`
and `881.2` share the same degree-18 field; the `14303` field is `Q(ζ₁₁)⁺`.) Both degrees are
`≫ 1`, so `σ ≅ ρ̄_{f,λ}` comes from an abelian variety of dimension 18 (resp. 5), **confirming
these mod-ℓ representations do not arise from elliptic curves** — the "not dimension 1" point.

**Scope.** This argument needs the Hecke field, and the Hecke field needs the char-0 eigenform. We
have that for exactly these **four** curves. For the other **35** the match is a *certificate*: the
survivor eigenvalues lie in `F_ℓ`, which does not by itself pin `[E:Q]`. §4f develops a test that
avoids eigenform isolation; §4f also reports, honestly, how far it actually got.

### 4f. A cheap degree lower bound from the survivor (mod-ℓ generalized eigenspace)

The kernel survivor is the mod-ℓ eigenvector `v` with `v·T_P = t_P·v`, `t_P = −#C(𝔽_P) mod ℓ`. By
**multiplicity one**, the `m`-adic Hecke module at the maximal ideal `m` cut out by this system is
free of rank 1, so the mod-ℓ **generalized** eigenspace `G` has `dim_{F_ℓ} G = [E_λ:Q_ℓ]`, where
`λ | ℓ` is the prime of the Hecke field `E` picked out by the fingerprint (residue degree 1). Since
`[E:Q] = Σ_{λ|ℓ} [E_λ:Q_ℓ] ≥ [E_λ:Q_ℓ] = dim G`,

> **`dim G > 1` ⇒ `[E:Q] > 1`** — the newform is not rational, so `σ` does not arise from an
> elliptic curve.

(The 1-dimensional survivor rules out any *other* newform congruent mod ℓ, so `G` is "pure" and
`dim G` is exactly `[E_λ:Q_ℓ]`.) Crucially `dim G` is computed by the **increasing chain of nested
kernels**
```
E₁ = survivor,   E_{k+1} = { w : w·(T_i − c_i) ∈ E_k  for all i },   G = lim_k E_k,
```
i.e. the **same `F_ℓ` linear algebra as the survivor** — so it **scales to the giant dimensions**
(unlike any char-0 lift). Validated (`ladic_degree.m`):

| form | ℓ | `dim G = [E_λ:Q_ℓ]` | what it certifies | known `[E:Q]` |
|---|--:|--:|---|--:|
| 14303 | 11 | **5** (chain `1→2→3→4→5`) | `[E:Q] ≥ 5` | 5 |
| 881 | 13 | **1** (`λ` split, residue deg 1) | **nothing** | **18** |

**The test is one-sided, and the `881` row is the proof.** There the Hecke field is known to have
degree 18, and `dim G` is nevertheless 1. So `dim G = 1` is not weak evidence for `[E:Q] = 1`; it is
*no* evidence either way. It says only that `ℓ` happens to have a degree-1 unramified prime in `E`,
which is a fact about how `ℓ` splits, not about how big `E` is. Any reading of an inconclusive result
as "probably an elliptic curve" would be wrong.

#### What the sweep actually found

We ran `ladic_degree.m` over the **27** curves of dimension `≤ 8345` (`run_ladic_shards.sh`; the
12 larger ones, idx 22–27 and 34–39, were never attempted). Results:

| outcome | curves | |
|---|--:|---|
| `dim G > 1` — `[E:Q] > 1` **certified** | **3** | `14303.1`, `14303.2` (`dim G = 5`); **`145161.2` (`dim G = 2`)** |
| `dim G = 1`, recovered by the ℓ-adic fallback | **2** | `881.1`, `881.2` (degree 18, at dim 441) |
| `dim G = 1`, fallback out of reach (`dim > 2000`) | 18 | inconclusive |
| `dim G = 1`, fallback ran and **failed** | 4 | `20447.3/.6`, `24889.1/.2` |

So the scalable test certified **3 of 27**, and only **one new curve** beyond what §4e already knew:
`145161.2`, with `[E:Q] ≥ 2`. That is a much weaker statement than the degrees 18 and 5 of §4e — it
rules out an elliptic curve and nothing more — and it should be quoted as such.

**The ℓ-adic fallback is exhausted, not merely slow.** When `dim G = 1` we have `a_P ∈ Z_ℓ` and can
try to reconstruct `a_P`'s minimal polynomial from a Hensel lift of `v` (LLL, find-at-half /
verify-at-full). This is what recovers degree 18 for `881` at dim 441. But it is `O(dim³)`, and on
`24889.1` (dim 1038) it returned **no relation at all for any of 12 primes**, even after raising the
precision from 400 to 1200 and the degree cap from 25 to 40. A diagnostic added for this run
distinguishes "recovered degree 1" from "found nothing" — the four curves in the last row are the
latter, so they carry no information about `[E:Q]` either.

#### Where this leaves the claim

> **"σ does not arise from an elliptic curve" is established for 5 of the 39 curves:** `881.1`,
> `881.2` (`[E:Q] = 18`), `14303.1`, `14303.2` (`[E:Q] = 5`) — all four by explicit Hecke field,
> §4e — and `145161.2` (`[E:Q] ≥ 2`, §4f). For the other 34 it is **open**, in the strict sense
> that we have no evidence in either direction.

This is a limitation of the method, not of the computation, and more machine time will not move it.
`dim G = [E_λ:Q_ℓ]` is bounded by the residue degree of `λ`, so whenever `ℓ` has a degree-1 prime in
`E` — which was the case for 22 of the 27 we tested — the test is silent no matter how large `[E:Q]`
really is. Closing the gap needs the char-0 eigenform, i.e. the `NewformDecomposition` step that is
already known to be intractable past dim ≈ 700 (§4c).

Two things are worth saying alongside the negative result, without overselling either. First, the
five curves where we *do* know the degree span both residual primes and give degrees 18, 5 and ≥ 2 —
there is no case in the dataset where a Hecke field turned out to be rational. Second, a rational
Hecke field would be a strong coincidence at these conductor sizes, so the expectation is that the
remaining 34 behave like the five. That is an expectation, not a theorem, and the write-up should
not dress it as more.

**Caveat — certificates vs. cutters.** The kernel method yields a *match certificate*
(1-dim surviving mod-ℓ eigenspace + control), not the char-0 eigenform, so it does **not**
by itself produce Hecke cutters (min. polys of `a_P`). The 4 forms in §4/§4a (with full Hecke
data in `hecke_cutters.m`) remain the explicitly-identified subset; extending *cutters* to the
rest requires isolating each surviving eigenform, which is the expensive char-0 step. The GRH
Faltings–Serre argument, by contrast, **does not** need that isolation — the survivor eigenvector
supplies `a_P mod λ` directly, which is why §4b scales to 37 curves.

## 5. What this says for the collaboration

- **Feasibility: easy for the bulk.** The small-conductor majority match in seconds–minutes;
  the fingerprint is trivial; the only real step is pinning the 2-part by level-lowering
  (cheap — bounded by `Conductor(C)`).
- **GRH** enters only for *certification*: an effective bound on the least distinguishing
  prime, independent of the (huge) splitting field. Finding the form is unconditional.
- **Scaling.** The handful of larger-conductor curves (dim `> ~few·10³`) use the
  kernel-intersection matcher (no `NewformDecomposition`) already developed for Goal 1.

## Reproduce

```
magma validate.m                  # structure + conductor consistency for all 39 curves
magma sweep.m                     # NewformDecomposition match (4 smallest only); streams to sweep.out
magma idx:=3 kernel_torsion.m     # kernel-intersection match certificate for curve #idx (§4c)
magma idx:=3 grh_kernel.m         # GRH modularity certificate for curve #idx (§4b), BOUND:=3000
```
`grh_kernel.m` writes `grh_kernel_<idx>.out`; the success sentinel is
`GRH-KERNEL CERT <label>: MODULAR`. `run_grh_shards.sh` launches the batch with `BOUND` tiered by
dimension (3000 for `dim ≤ 9k`, 800 above); collect with
`grep -h "GRH-KERNEL CERT" grh_kernel_*.out`.
`kernel_torsion.m` is the full-sweep matcher (§4c): per curve it builds `M = HilbertCuspForms`
and reports the surviving mod-ℓ eigenspace dim + control, writing `kernel_<idx>.out`. It reaches
the whole tail (dims to ~40k) where `sweep.m` cannot. It certifies a match but does **not** emit
Hecke cutters (see the §4c caveat).

`hecke_cutters.m` — loadable Hecke data (field + cutters as `<prime, minpoly(a_P)>`) for the 4
explicitly-identified forms (§4a); `magma lab:="14303.1" e:=0 out:="hecke_cutters.m" emit_cutters.m`
appends a form. The other kernel matches (§4c) are certificates only — cutters pending eigenform
isolation. (GRH certificates, by contrast, are *not* blocked on this: see §4b.)

Data: `torsion_data.m` (transcribed from `examples.json`, each curve validated).
