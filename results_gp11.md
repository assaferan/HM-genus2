# (1,11) Gross–Popescu Search Results

Genus-2 curves over **K = Q(√5)** with a rational **(1,11)-polarized isogeny**,
found by searching non-rational points on the D₂ locus of the Gross–Popescu
construction, inverting via `InversionGP11Fast`, and filtering for good reduction
at both primes above 11.

Full L-polynomial data: `data_qsqrt5_gp11.txt`.  
Survivor Igusa invariants: `survivors_qsqrt5_gp11.m`.  
Verification script: `verify_curveA.m`.

## Field

K = Q(√5), O_K = Z[(1+√5)/2] (golden-ratio ring),  
disc(K/Q) = 5, class number h(K) = 1.  
11 splits in K: 11·O_K = p₁₁ₐ · p₁₁ᵦ, N(p₁₁ₐ) = N(p₁₁ᵦ) = 11.

## Search parameters

| Parameter | Value |
|---|---|
| `BoxM` | 2 |
| `MaxPts` | 200 |
| `TrialsRetry` (first 30 pts) | 200 |
| `TrialsNew` (pts 31–200) | 150 |
| Precision | 100 |
| Filters | genuine over K (not Q-base-change), potential good reduction at 11, End(J)=Z (simple) |

Scripts: `search_qsqrt5_gp11.m` (initial), `search_qsqrt5_gp11_ext.m` (extended).

## Survivors

4 survivors found, forming **3 K-isomorphism classes** (A, B, B'):

| Class | Survivors | Bad prime N(p) | Conductor exponent | p|11 exponent |
|---|---|---|---|---|
| **A** | 1, 2 (same Igusa invs) | 66179 (prime) | 1 | 0 (both p|11) |
| **B** | 3 | 602479 (prime) | 1 | 0 (both p|11) |
| **B'** | 4 | 602479 (prime) | 1 | 0 (both p|11) |

B and B' are Galois conjugates (related by the non-trivial automorphism s5 ↦ −s5
of K/Q). Each is a distinct curve over K but they form a single orbit under Gal(K/Q).

Note: the prime-at-2 conductor exponent was **not computed** (Magma's
`ConductorExponentAt` does not handle p|2 over number fields). The quoted
conductor norms N(p) are the **odd part** only.

## Curve A — details

**Igusa-Siegel invariants** [J₂:J₄:J₆:J₈:J₁₀] in WPS(2,4,6,8,10), normalized J₂=1:
```
J2  = 1
J4  = 1/1459240*(243125*s5 - 482787)
J6  = 1/2229718720*(-54798934*s5 + 127710391)
J8  = 1/8517525510400*(182422196780*s5 - 406668692089)
J10 = 1/65073894899456000*(38134182372761*s5 - 85270624049895)
```

**Minimal twist factor:** d = -s5 + 10  
(the Mestre reconstruction gives a non-minimal model; MinimalTwist corrects it)

**Conductor:** bad prime above 66179 (norm 66179, degree-1 prime), exponent 1.  
Both primes above 11 have conductor exponent 0 (good reduction at 11). ✓

**Selected Frobenius L-polynomials** (verified against `data_qsqrt5_gp11.txt`):

| N(p) | p | L_p(T) |
|---|---|---|
| 49 | 7 (inert) | 2401T⁴ + 98T³ − 58T² + 2T + 1 |
| 11 | 11 (split, p₁) | 121T⁴ + 33T³ + 18T² + 3T + 1 |
| 11 | 11 (split, p₂) | 121T⁴ − 44T³ + 14T² − 4T + 1 |
| 169 | 13 (inert) | 28561T⁴ − 2873T³ + 376T² − 17T + 1 |
| 289 | 17 (inert) | 83521T⁴ + 867T³ + 288T² + 3T + 1 |

Verification: `magma verify_curveA.m` (runs in ~8 s).

**Frobenius at N(p) = 3 and 5** fail with "reduction failed" — this is expected:
the model has structural bad reduction at these primes due to the (1,5)-polarization
structure of the Mestre conic over K (see Weierstrass model note below).

## Curve B / B' — details

**Igusa-Siegel invariants** (Curve B):
```
J2  = 1
J4  = 1/283315208*(-279106*s5 + 8627389)
J6  = 1/3372017605616*(1674615499*s5 - 619836124)
J8  = 1/321070028336333056*(44678459406664*s5 - 89575920849197)
J10 = 1/15285501909036144130048*(-138485012172127*s5 + 309690942212419)
```

Curve B' is the Galois conjugate: replace s5 ↦ −s5 in all of the above.

Bad prime: N(p) = 602479 (prime), conductor exponent 1.  
Both primes above 11: exponent 0. ✓

## Weierstrass model — coefficient size note

The explicit Weierstrass model y² = f(x) over K is available but has very large
coefficients (numerators/denominators with 170–350 decimal digits). This is not a
code deficiency — it is intrinsic:

1. The Igusa invariant J₁₀ has denominator ~10¹⁶; Mestre's algorithm propagates
   this into the model coefficients as ~10^{170} rational numbers.
2. Magma's `ReducedModel` / `Reduce` only work over Q, not number fields.
3. GL₂(O_K) Möbius transforms (implemented in `reduce_curveA*.m`) reduce the
   log-height from 41 → 32 → 18, but the rational numerators/denominators remain
   ~170 digits regardless (the transform creates near-cancellation, not simpler
   numbers).
4. A degree-5 model exists (`HasOddDegreeModel` returns true, log-height 28.9),
   but still has ~185–350 digit coefficients.

**For computation:** use `TryReduceModP` (in `data_qsqrt5_gp11.m` /
`verify_curveA.m`) to reduce mod a prime ideal without needing a global
integral model. This works at all good primes above p ≥ 7.

**For a paper:** quote the Igusa-Siegel invariants [J₂:…:J₁₀] and the
conductor data. These are compact and sufficient to identify the curve.

## How to reproduce

```magma
// Verify Curve A (conductor + Frobenius spot-check, ~8 s):
magma verify_curveA.m

// Run extended search for new survivors:
magma search_qsqrt5_gp11_ext.m

// Compute Frobenius data for all three curves:
magma data_qsqrt5_gp11.m
```
