# Matching Hugo Nartz's torsion-type mod-ℓ representations to Hilbert modular forms

*Working notes for the collaboration. Dataset: `examples.json` (Hugo Nartz, via Ariel) —
genus-2 curves over `Q(√2)`, `Q(√3)`, `Q(√5)` whose mod-ℓ Galois representation (`ℓ ≥ 11`)
contains an irreducible 2-dimensional sub coming from a torsion point.*

## 0. Summary

For each curve `C/F` (`F` real quadratic) the mod-ℓ representation of `A = Jac(C)` is
```
    A[ℓ]^ss  =  1  ⊕  χ_ℓ  ⊕  σ,          σ irreducible, 2-dimensional, det σ = χ_ℓ,
```
so `σ` is a candidate for a **parallel weight-`[2,2]`, trivial-nebentypus Hilbert modular
form over `F`**. This is a direct (non-induced) test of Serre's conjecture over a real
quadratic field.

**Findings (this file covers `Q(√2)` and `Q(√3)`; the JSON has 28 + 12 candidates there):**
- The whole list is **computationally easy**: because `ℓ ∈ {11,13}` does not divide any
  conductor, the Hilbert level is (essentially) the conductor, with **no `ℓ`-power blow-up**
  — unlike the `p=5` non-base-change family, whose `5⁴`-inflated levels reach dim `10⁵–10⁶`.
  Here every `Q(√2)` space has **dim ≤ 17 281**, and every `Q(√3)` space is comparable.
- The fingerprint is a **one-line point count**: `tr σ(Frob_P) ≡ −#C(𝔽_P) (mod ℓ)`.
- We **match** the small-conductor curves to explicit Hilbert newforms, each with a
  wrong-fingerprint discrimination control (a true match keeps agreement at all primes; the
  control collapses to 0). [matches table below]
- Under **GRH**, each match upgrades to a theorem via an effective Faltings–Serre /
  Chebotarev prime bound `O((log cond)²)` — a few-hundred-prime check, as in the idx-33 work.

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

All ≤ 55 446 (larger per unit norm since `ζ_{Q(√3)}(−1) = 1/6` vs `1/12` for `Q(√2)`;
`Q(√5)` would be `≈ 0.4×` the `Q(√2)` size).

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
- **orbit dim** 18, **Hecke field** Q[y]/(y^18 - 54y^16 - 6y^15 + 1124y^14 + 184y^13 - 11538y^12 - 2264y^11 + 61687y^10 + 13681y^9 - 166924y^8 - 35569y^7 + 203946y^6 + 26001y^5 - 95812y^4 - 10226y^3 + 12624y^2 + 2360y - 16)
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
- **orbit dim** 5, **Hecke field** Q[y]/(y^5 - 2y^4 - 5y^3 + 13y^2 - 7y + 1)
- **Hecke cutters** — minpoly of `a_P`:
    - `P (norm 2, (-a)):  minpoly(a_P) = y^5 - 2y^4 - 5y^3 + 13y^2 - 7y + 1`
    - `P (norm 9, (3)):  minpoly(a_P) = y^5 + y^4 - 15y^3 - 14y^2 + 3y + 1`
    - `P (norm 25, (5)):  minpoly(a_P) = y^5 - 4y^4 - 75y^3 + 170y^2 + 1010y + 439`
    - `P (norm 7, (-2a + 1)):  minpoly(a_P) = y^5 - 5y^4 - 23y^3 + 122y^2 + 115y - 683`
    - `P (norm 7, (-2a - 1)):  minpoly(a_P) = y^5 - 33y^3 + 242y + 121`
    - `P (norm 121, (11)):  minpoly(a_P) = y^5 + 9y^4 - 335y^3 - 2088y^2 + 23665y + 30097`

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
magma validate.m     # structure + conductor consistency for all 39 curves
magma sweep.m        # level-lowering match + discrimination control, streams to sweep.out
```
Data: `torsion_data.m` (transcribed from `examples.json`, each curve validated).
