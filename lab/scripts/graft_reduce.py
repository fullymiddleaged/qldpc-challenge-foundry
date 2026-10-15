# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Shrink a CSS code by row-space grafts (qec_search.graft, @dorakingx's move) at fixed k, weight cap, radius cap and
distance, checking each non-free graft with upstream's gf2_fast RIS: a --screen-trials pass, then --confirm-rounds
passes of --confirm-trials on fresh seeds, all required to find nothing lighter than --d. Those passes are a filter;
the evidence for the reduced code is a fresh ladder (scripts/ris_ladder.py) afterwards. A 2k + 20k filter let a
distance-losing graft through at n = 392 (15 Oct 2026: the result read d <= 16 at 200k), hence the defaults.

Run in WSL from lab/ (gf2_fast is built there, see scripts/ris_ladder.py):
    python3 scripts/graft_reduce.py results/hunt/a637_scale/cand/a637_10x10.npz --d 15 --gf2 ~/build/gf2_fast
writes <npz stem>_graft.npz (hx, hz, coords, idx = source columns kept) and <stem>_graft.log beside it, or
--out <path>.npz and <path>.log.
"""
from __future__ import annotations

import argparse
import os
import pathlib
import sys
import time
import zlib

import numpy as np

LAB = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB))
from qec_search import graft  # noqa: E402


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("npz")
    p.add_argument("--d", type=int, required=True, help="distance the reduced code must keep")
    p.add_argument("--gf2", required=True, help="directory holding the built gf2_fast module")
    p.add_argument("--w", type=int, default=8)
    p.add_argument("--radius", type=float, default=4.0)
    p.add_argument("--terms", type=int, default=3, help="S is a sum of at most this many checks")
    p.add_argument("--screen-trials", type=int, default=20000)
    p.add_argument("--confirm-trials", type=int, default=200000)
    p.add_argument("--confirm-rounds", type=int, default=2)
    p.add_argument("--final-trials", type=int, default=1000000,
                   help="at the fixpoint, two passes of this many trials; on failure roll back to the last state that passes")
    p.add_argument("--out", type=pathlib.Path, help="output .npz (log beside it); default <npz stem>_graft.npz")
    p.add_argument("--threads", type=int, default=os.cpu_count())
    a = p.parse_args(argv)
    sys.path.insert(0, os.path.expanduser(a.gf2))
    import gf2_fast

    src = pathlib.Path(a.npz)
    z = np.load(src)
    hx, hz, coords = z["hx"].astype(np.uint8), z["hz"].astype(np.uint8), z["coords"].astype(float)
    out = a.out or src.with_name(src.stem + "_graft.npz")
    logf = out.with_suffix(".log")
    calls = [0]

    def passes(nx, nz, budgets) -> bool:
        # hashed seeds: consecutive small seeds (1, 2, 3, ...) explored poorly and passed a d-losing graft at n = 512
        for t in budgets:
            calls[0] += 1
            s = zlib.crc32(f"{src.name}|graft|{calls[0]}".encode()) & 0x7FFFFFFF
            if gf2_fast.distance_rand_parallel(nx, nz, trials=t, seed=s, threads=a.threads) < a.d:
                return False
        return True

    def dist_ok(nx, nz) -> bool:
        return passes(nx, nz, (a.screen_trials, *[a.confirm_trials] * a.confirm_rounds))

    t0 = time.time()
    with logf.open("w", encoding="utf-8") as f:
        def log(msg: str) -> None:
            line = f"[{time.time() - t0:7.0f}s] {msg}"
            print(line, flush=True)
            f.write(line + "\n")
            f.flush()

        log(f"source {src.name}: n={hx.shape[1]} k={graft.k_of(hx, hz)} d>={a.d} w<={a.w} radius<={a.radius}")
        states = [(hx, hz, coords, np.arange(hx.shape[1]))]
        graft.reduce(hx, hz, coords, a.w, a.radius, dist_ok, a.terms, log, on_step=lambda *st: states.append(st))
        log(f"fixpoint after {len(states) - 1} grafts; final check {a.final_trials} trials x 2")
        best = graft.last_good(len(states), lambda i: passes(*states[i][:2], (a.final_trials,) * 2))
        if best < len(states) - 1:
            log(f"final check failed: rolled back to graft {best} of {len(states) - 1}")
        rx, rz, rc, idx = states[best]
        np.savez_compressed(out, hx=rx, hz=rz, coords=rc, idx=idx)
        log(f"fixpoint: n={rx.shape[1]} k={graft.k_of(rx, rz)}; wrote {out.name}")


if __name__ == "__main__":
    main()
