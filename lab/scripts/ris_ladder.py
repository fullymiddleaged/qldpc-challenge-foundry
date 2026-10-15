# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Fresh-seed RIS ladder on CSS codes with upstream's bit-packed C++ engine (verify/gf2_fast.cpp).

Each rung runs gf2_fast.distance_rand_witness with its own seed and appends one JSON line (code, trials, seed,
d_ub, side, witness, seconds) to --log. Rungs already in the log are skipped, so a killed run resumes where it
stopped. A rung that finds a logical below an earlier rung's weight lowers the code's d_ub; later rungs still run.
The result is an upper bound, never a proof.

gf2_fast does not build on this Windows box (no C++ compiler); build it in WSL once:
    python3 verify/setup_gf2_fast.py build_ext --build-lib ~/build/gf2_fast --build-temp /tmp/gf2tmp
then run from lab/ in WSL:
    python3 scripts/ris_ladder.py --gf2 ~/build/gf2_fast --log results/hunt/a637_scale/ladder.jsonl \\
        --ladder 20000,200000,1000000,5000000 results/hunt/a637_scale/cand/*.npz
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time
import zlib

import numpy as np


def done_rungs(log: pathlib.Path) -> set[tuple[str, int]]:
    if not log.exists():
        return set()
    out = set()
    for line in log.read_text(encoding="utf-8").splitlines():
        if line.strip():
            r = json.loads(line)
            out.add((r["npz"], r["trials"]))
    return out


def rung_seed(name: str, trials: int) -> int:
    """Distinct per code and rung, stable across restarts."""
    return zlib.crc32(f"{name}|{trials}".encode()) & 0x7FFFFFFF


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("npz", nargs="+")
    p.add_argument("--gf2", required=True, help="directory holding the built gf2_fast module")
    p.add_argument("--log", required=True)
    p.add_argument("--ladder", default="20000,200000,1000000,5000000")
    p.add_argument("--threads", type=int, default=os.cpu_count())
    a = p.parse_args(argv)
    sys.path.insert(0, os.path.expanduser(a.gf2))
    import gf2_fast

    log = pathlib.Path(a.log)
    done = done_rungs(log)
    ladder = [int(t) for t in a.ladder.split(",")]
    for path in a.npz:
        name = pathlib.Path(path).as_posix()
        rel = name.split(log.parent.as_posix() + "/", 1)[-1]
        z = np.load(path)
        hx, hz = z["hx"].astype(np.uint8), z["hz"].astype(np.uint8)
        for t in ladder:
            if (rel, t) in done:
                continue
            seed = rung_seed(rel, t)
            t0 = time.time()
            w, side, supp = gf2_fast.distance_rand_witness(hx, hz, trials=t, seed=seed, threads=a.threads)
            rec = {"npz": rel, "n": int(hx.shape[1]), "trials": t, "seed": seed, "d_ub": w, "side": side,
                   "witness": [int(q) for q in supp], "secs": round(time.time() - t0, 1)}
            with log.open("a", encoding="utf-8") as f:
                f.write(json.dumps(rec) + "\n")
            print(json.dumps({k: v for k, v in rec.items() if k != "witness"}), flush=True)


if __name__ == "__main__":
    main()
