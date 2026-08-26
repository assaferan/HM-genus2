#!/usr/bin/env python3
"""Check torsion_data.m against the source dataset (Hugo Nartz's examples.json).

Everything in this repo rests on torsion_data.m being a faithful transcription of the
supplied curves: a wrong coefficient would not announce itself -- the pipeline would
happily match the *wrong* curve to a Hilbert newform and certify it. So verify rather
than trust.

Compares, per curve and in order: base field d, residual prime l, ideal label,
prime-to-2 conductor norm, and both coefficient lists (f and h, each [u,v] = u + v*a).

Usage:  ./verify_transcription.py [path/to/examples.json]

examples.json is not committed (it is collaborators' data, and torsion_data.m already
encodes all of it). Supply it locally to run the check; run_tests.sh skips this gate
when the file is absent.

Sentinel on success: "TRANSCRIPTION: ALL PASS (n/n curves)"
"""
import ast
import json
import os
import re
import sys

FIELD_D = {"Q(sqrt(2))": 2, "Q(sqrt(3))": 3, "Q(sqrt(5))": 5}


def load_source(path):
    with open(path) as fh:
        blob = json.load(fh)
    out = []
    for fld in blob["fields"]:
        name = fld["field"]
        if name not in FIELD_D:
            raise SystemExit(f"unknown field in dataset: {name}")
        for cand in fld["candidates"]:
            cond = cand["abelian_variety_conductor_prime_to_2"]
            out.append(
                dict(d=FIELD_D[name], l=cand["residual_prime"], lab=cond["ideal_label"],
                     N=cond["norm"], fc=cand["equation"][0], hc=cand["equation"][1])
            )
    return out


def _read_list(s, i):
    """Read one balanced [...] literal starting at s[i]; return (value, index_after)."""
    while s[i] in " ,\n\t":
        i += 1
    if s[i] != "[":
        raise ValueError(f"expected '[' at {i}: {s[i:i+30]!r}")
    depth = 0
    for j in range(i, len(s)):
        if s[j] == "[":
            depth += 1
        elif s[j] == "]":
            depth -= 1
            if depth == 0:
                return ast.literal_eval(s[i:j + 1]), j + 1
    raise ValueError("unbalanced brackets")


def load_transcription(path):
    txt = open(path).read()
    out = []
    for m in re.finditer(r"\bd([235])\(", txt):
        d = int(m.group(1))
        i = m.end()
        head = re.match(r'\s*(\d+)\s*,\s*"([^"]+)"\s*,', txt[i:])
        if not head:
            continue
        l, lab = int(head.group(1)), head.group(2)
        i += head.end()
        fc, i = _read_list(txt, i)
        hc, i = _read_list(txt, i)
        tail = re.match(r"\s*,\s*(\d+)", txt[i:])
        if not tail:
            raise SystemExit(f"could not parse conductor norm for {lab}")
        out.append(dict(d=d, l=l, lab=lab, N=int(tail.group(1)), fc=fc, hc=hc))
    return out


def main():
    src_path = sys.argv[1] if len(sys.argv) > 1 else "examples.json"
    here = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(src_path):
        print(f"TRANSCRIPTION: SKIP (no {src_path}; supply the dataset to run this check)")
        return 0

    src = load_source(src_path)
    got = load_transcription(os.path.join(here, "torsion_data.m"))

    print(f"  source dataset:  {len(src)} curves")
    print(f"  torsion_data.m:  {len(got)} curves")
    if len(src) != len(got):
        print(f"TRANSCRIPTION: FAIL (count mismatch)")
        return 1

    bad = 0
    for k, (a, b) in enumerate(zip(src, got), 1):
        for f in ("d", "l", "lab", "N", "fc", "hc"):
            if a[f] != b[f]:
                print(f"  FAIL idx {k} ({a['lab']}): {f}\n       source: {a[f]}\n       repo:   {b[f]}")
                bad += 1

    by_field = {}
    for c in src:
        by_field[c["d"]] = by_field.get(c["d"], 0) + 1
    print("  fields: " + ", ".join(f"Q(sqrt {d}): {n}" for d, n in sorted(by_field.items())))

    if bad:
        print(f"TRANSCRIPTION: FAIL ({bad} mismatched fields)")
        return 1
    print(f"TRANSCRIPTION: ALL PASS ({len(src)}/{len(src)} curves)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
