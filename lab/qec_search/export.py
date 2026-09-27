# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Export the best distinct codes from a run as .npz files for the
Unitary Foundation qLDPC Challenge submission tool (`./qldpc submit x.npz`).

    python -m qec_search.export runs/run1 --top 6 --dest submissions/

Picks codes that are not dominated by any other code in the run on
(n lower, k higher, d higher), after dropping split codes and relabelled
duplicates. Writes <dest>/<name>.npz with arrays hx, hz, and a
manifest.json describing each code's construction.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .bbcode import BBGenome, build_bb
from .known_codes import KNOWN_SIGS


def pareto(rows):
    front = []
    for r in rows:
        dominated = any(o is not r and o["n"] <= r["n"] and o["k"] >= r["k"] and o["d_ub"] >= r["d_ub"]
                        and (o["n"], o["k"], o["d_ub"]) != (r["n"], r["k"], r["d_ub"]) for o in rows)
        if not dominated:
            front.append(r)
    return front


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--dest", default="submissions")
    ap.add_argument("--layout", action="store_true",
                    help="also search for a 2D bilayer layout (max check diameter <= 7) for each code")
    args = ap.parse_args(argv)

    rows = {}
    for d in args.runs:
        p = Path(d) / "results.jsonl"
        if not p.exists():
            continue
        for line in p.open():
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get("k") and r.get("d_ub") and not r.get("split") and not r.get("dup_of"):
                rows.setdefault(r["key"], r)

    # re-check split / duplicates for records made by older versions
    by_sig = {}
    for r in sorted(rows.values(), key=lambda r: -r.get("kd2n", 0)):
        code = build_bb(BBGenome.from_json(r["genome"]))
        if code.components() > 1:
            continue
        sig = code.signature()
        if sig in KNOWN_SIGS:  # an IBM code, possibly relabelled
            continue
        if sig not in by_sig:
            by_sig[sig] = r
    distinct = list(by_sig.values())
    # IBM's published codes are not ours to submit; and one code per
    # [[n,k,d]] is enough, since the board ranks on those numbers only
    # (keep the one that decoded best in our Monte Carlo)
    def ler(r):
        cc = r.get("cc") or []
        return cc[-1]["ler"] if cc else 1.0
    per_params = {}
    for r in distinct:
        if r.get("known"):
            continue
        key = (r["n"], r["k"], r["d_ub"])
        if key not in per_params or ler(r) < ler(per_params[key]):
            per_params[key] = r
    front = sorted(pareto(list(per_params.values())), key=lambda r: -r["kd2n"])[: args.top]

    dest = Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for r in front:
        g = BBGenome.from_json(r["genome"])
        code = build_bb(g)
        name = f"bb_{code.n}_{code.k}_{r['d_ub']}_{g.l}x{g.m}"
        arrays = {"hx": code.hx.astype(np.uint8), "hz": code.hz.astype(np.uint8)}
        layers, radius = 1, None
        if args.layout:
            from .layout import find_layout, validate, check_supports
            for pad in (2, 5):
                coords, _ = find_layout(code, g, target=7.0, steps=2_000_000, restarts=2, pad=pad)
                v = validate(coords, check_supports(code))
                if v["radius"] <= 7.0 and v["ok_sites"]:
                    arrays["coords"], layers, radius = coords, 2, round(v["radius"], 3)
                    break
        np.savez(dest / f"{name}.npz", **arrays)
        manifest[name] = {
            "layers": layers,
            "radius": radius,
            "params": f"[[{code.n},{code.k},{r['d_ub']}]]",
            "kd2n": r["kd2n"],
            "construction": f"Bivariate bicycle code on Z_{g.l} x Z_{g.m}: {g.pretty().split(', ', 2)[2]}. "
                            f"Found by evolutionary search (qec-search)"
                            + (f"; bilayer layout (max check diameter {radius}) by simulated annealing." if radius else "."),
            "genome": r["genome"],
        }
        print(f"{name:28s} [[{code.n},{code.k},{r['d_ub']}]]  k*d^2/n={r['kd2n']:.2f}  "
              + (f"bilayer layout radius {radius}  " if radius else ("no layout found  " if args.layout else ""))
              + g.pretty(), flush=True)
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"\n{len(front)} codes written to {dest} (Pareto front of {len(distinct)} distinct codes)")


if __name__ == "__main__":
    main()
