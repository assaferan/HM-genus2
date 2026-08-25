# Magma #110 workaround — definite Hilbert modular forms basis inverse

Three of the 39 torsion-type curves (the largest `Q(√3)` spaces) cannot be
matched with stock Magma because building their Hilbert cusp space crashes:

| idx | label      | ℓ  | level  | dim ≈ |
|-----|------------|----|--------|-------|
| 37  | 569399.3   | 13 | 569399 | 47168 |
| 38  | 785473.12  | 11 | 785473 | 55446 |
| 39  | 785473.5   | 11 | 785473 | 55446 |

## The bug

On the first `HeckeOperator` call, `BasisMatrixDefinite` (in
`package/Geometry/ModFrmHil/definite.m`) throws:

```
BasisMatrixDefinite(M) → definite.m line 1060:
    Binv := Transpose(Solution(Transpose(B), IdentityMatrix(BaseRing(B), Nrows(B))));
Runtime error in 'Solution': No solution exists
```

Deterministic, and **not** a size limit — dim 39418 (472993) succeeds.

## Root cause

`B` (`basis_matrix_big`) is assembled fine from the ideal-class direct factors.
The crash is one step later: Magma computes a right inverse `Binv` of `B` via
`Solution(Bᵀ, I)`, which requires `B` to have full **row** rank. For these two
levels the assembly yields **linearly dependent rows**, so the solve has no
solution. (A genuine Magma bug, filed as Magma issue #110.)

## Why the workaround is safe for parallel weight 2

`basis_matrix_big_inv` (the `Binv` computed at line 1060) is **read in exactly
one place** — `definite.m:1073`, inside the *non*-weight-2 branch. For weight
`[2,2]`:

* `RemoveEisenstein` rebuilds `basis_matrix` / `basis_matrix_inv` from the
  Eisenstein indicator vectors and the inner product — never touching `Binv`;
* the big Hecke matrix uses only `Ncols(basis_matrix_big)`, not the inverse.

So for parallel weight 2 the crashing `Binv` is **write-only / vestigial**. The
patch (`definite.m.patch`) wraps the solve in `try/catch`: on failure it skips
`Binv` for weight 2, and **re-raises for any other weight** (where `Binv` is
consumed). On the success path it is byte-for-byte the original computation.

## Validation

Patched vs stock Magma, kernel-intersection matcher (`kernel_torsion.m`):

| curve            | field   | result                                | vs stock        |
|------------------|---------|---------------------------------------|-----------------|
| idx 3 · 14303.1  | Q(√2)   | MATCH e=0, survivor=1/control=0       | **identical**   |
| idx 1 · 881.1    | Q(√2)   | reject e=0,1,2 → MATCH e=3             | matches known   |
| idx 28 · 4057.1  | Q(√3)   | reject e=0,1 → MATCH e=2               | **identical**   |

The patch is inert wherever line 1060 succeeds; it only changes the behaviour
of the previously-crashing giants.

## Deploy (no system files touched)

```
./deploy_patch.sh <magma-install-dir> <dest-dir>
# e.g. on lovelace:
./deploy_patch.sh /opt/magma/magma-2.29-9 /scratch/$USER/magma-patched
```

This copies a Magma install you already have a licence for and applies our diff
to the copy. Run `<dest-dir>/magma` for the affected computations. We never
redistribute Magma source — only `definite.m.patch` (our ~25 lines) is tracked
here; the full files are gitignored.
