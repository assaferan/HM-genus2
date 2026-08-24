# Banked evidence — matches, GRH certificates, and the ℓ-adic degree sweep

Raw output of the runs behind `results_torsion_hmf.md`. These are the product of several weeks
of compute on machines we do not control, and the individual runs are not reproducible on demand
(the largest cost ~38 h each and >400 GB). They are banked here so the results survive
independently of any scratch directory, and so they can be re-checked without re-running anything.

Run `../verify_evidence.sh` to validate the whole set against `torsion_data.m`
(pure text checking, no Magma, <1 s). It is part of `run_tests.sh`.

| directory | files | what |
|---|--:|---|
| `matches/` | 39 | `kernel_torsion.m` output: the kernel-intersection match per curve (§4c/§4d) |
| `certificates/` | 36 | `grh_kernel.m` output: the GRH modularity certificate per curve (§4b) |
| `ladic/` | 29 | `ladic_degree.m` output: the generalized-eigenspace degree sweep (§4f) |

**Status: 39/39 matched, 36/39 certified.** The three outstanding certificates are the §4d giants
(idx 37/38/39); those runs are in progress and this directory should be refreshed when they land.

## What each file proves

- `matches/kernel_<idx>.out` — ends in `KERNEL SHARD <idx> DONE` and carries a line
  `<label> ... survivor=1 control=0 MATCH (e=<e>)`. `survivor=1` is the 1-dimensional surviving
  mod-ℓ eigenspace; `control=0` is the wrong-fingerprint control collapsing, which is what makes
  the match meaningful rather than a coincidence.
- `certificates/grh_kernel_<idx>.out` — ends in `GRH-KERNEL CERT <label>: MODULAR` with
  `disagree=0` and `irred`, plus the count of good primes verified and the `BOUND` used.
- `ladic/*.log` — per-curve `dim G`, and (where the fallback ran) per-prime recovered degrees.
  Note that most of these are **inconclusive**, by the nature of the test; see §4f.

## Provenance, and a trap worth knowing about

The evidence was produced on two hosts that **share no filesystem**, and the authoritative copy of
a given file is not always where you would guess:

| idx | produced on | Magma |
|---|---|---|
| matches 1–36, certificates 1–36 | `lovelace` | stock V2.29-9 |
| matches 37–39 | `lovelace` (37) / `hensel` (38, 39) | patched V2.29-9 (Magma #110, see `../magma110_patch/`) |
| ladic sweep | `lovelace` | stock V2.29-9 |

**The trap:** `kernel_37.out`, `kernel_38.out` and `kernel_39.out` *also* exist in the repo
directory on lovelace, but those are **stale headers** — the first attempts at the giants were
OOM-killed before writing any result, and `kernel_torsion.m` truncates its output file at startup.
They look superficially fine (right name, plausible size) but contain no `MATCH` line. The real
idx-37 record lived only in a scratch workspace (`patch_validate/`) and was nearly lost.

`verify_evidence.sh` exists to make that class of error impossible to commit: it requires the
completion sentinel (which embeds the shard index, so transposed conjugates are caught too) and
checks the label inside each file against `torsion_data.m` by index. Both failure modes are
covered by negative tests.

## Regenerating

```
magma idx:=<i> kernel_torsion.m      # match           -> kernel_<i>.out
magma idx:=<i> BOUND:=800 grh_kernel.m   # certificate -> grh_kernel_<i>.out
magma idx:=<i> ladic_degree.m        # degree bound    -> stdout
```

The giants (idx 37/38/39) additionally need the `magma110_patch/` build on any Magma before
V2.29-10, and must be run **one at a time** — see `run_giants_sequential.sh` and
`run_grh_giants.sh`, whose headers record the measured costs.
