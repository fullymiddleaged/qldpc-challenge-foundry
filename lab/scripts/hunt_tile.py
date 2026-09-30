# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Hunt single-layer (interaction radius <= 4.0) tile codes (qec_search.tile, arXiv:2504.09171).

The paper optimised k d^2 / n with no locality cap and searched B = 3 exhaustively and B = 4 at random; the board's
tile entries (e-eight) are bilayer weight-8 at large n. Nobody has targeted the single-layer cell, whose frontier is
weak (best k d^2 / n about 5 on 29 Sep), and about half of random B = 3 weight-6 tiles already fit radius 4 in the
rotated layout (tile.layout, one qubit per site). The layout comes with the tile, so there is no annealing stage.

  screen  every tile of weight w in a B x B box (exhaustive for B <= 3, --samples at random above) whose X- and
          Z-tile radius is <= --cap; keep bulk sizes l x m (n <= --n-max) with k >= --kmin whose (n, k, d_ub, w) no
          entry on the board's single-layer weight-w frontier dominates. Every tile screened goes to screened.txt.
  prove   exact d via qec_search.distqldpc for the best --max by k d_ub^2 / n, with a board-twin check; writes an
          .npz (hx, hz, coords) per code.

    python scripts/hunt_tile.py screen --out results/hunt/tile1 --B 3 --w 6
    python scripts/hunt_tile.py prove --out results/hunt/tile1 --max 20
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import os
import pathlib
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

LAB = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB))
from qec_search import gf2, tile  # noqa: E402
from qec_search.bbcode import CSSCode  # noqa: E402
from qec_search import frontier as fr  # noqa: E402

_spec = importlib.util.spec_from_file_location("hunt_bilayer", LAB / "scripts" / "hunt_bilayer.py")
hb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hb)                      # board_twin, read_jsonl


def tile_radius(x_tile, B: int) -> float:
    """Largest check diameter a tile can give in tile.layout: its own X or Z support (boundary checks are subsets)."""
    best = 0.0
    for t in (x_tile, tile.z_tile(x_tile, B)):
        c = tile.layout(sorted(t)).astype(float)
        if len(c) > 1:
            d = np.sqrt(((c[:, None, :] - c[None, :, :]) ** 2).sum(-1))
            best = max(best, float(d.max()))
    return best


def tile_key(x_tile, B: int) -> str:
    return f"B{B}|" + ",".join(f"{k}{x}{y}" for k, x, y in sorted(x_tile))


def key_hash(key: str) -> str:
    return hashlib.sha1(key.encode()).hexdigest()[:12]


def all_tiles(B: int, w: int, samples: int, seed: int):
    edges = [(k, x, y) for k in "hv" for x in range(B) for y in range(B)]
    if B <= 3:
        yield from (list(t) for t in itertools.combinations(edges, w))
    else:
        rng = random.Random(seed)
        for _ in range(samples):
            yield rng.sample(edges, w)


def bulk_sizes(B: int, n_max: int, lmin: int = 3, aspect: float = 2.0):
    """l <= m <= aspect * l with n <= n_max. Long strips are capped by their short side and cost the most to screen
    (on 30 Sep a full sweep was 183 sizes and 260 s per tile)."""
    for l in range(lmin, 40):
        for m in range(l, int(aspect * l) + 1):
            if 2 * (l + B - 1) * (m + B - 1) <= n_max:
                yield l, m


def probe_size(B: int, n_probe: int) -> tuple[int, int]:
    """The square bulk whose n is closest to n_probe."""
    return min(((s, s) for s in range(3, 40)), key=lambda lm: abs(2 * (lm[0] + B - 1) ** 2 - n_probe))


def screen_one(x_tile, B: int, w: int, frontier, kmin: int, n_max: int, trials: int, cap: float,
               n_probe: int = 200, aspect: float = 2.0) -> list[dict]:
    """Radius, then connectivity and k at 6 x 6, then one probe size, then the sweep. A tile whose probe-size code is
    beaten is dropped without a sweep: a heuristic (a tile good only at other sizes is lost) that cuts the work about
    100x; loosen it with --n-probe 0 (no probe)."""
    if tile_radius(x_tile, B) > cap + 1e-9:
        return []
    hx, hz, _ = tile.build(x_tile, B, 6, 6)                    # cheap first look at one size
    code = CSSCode(hx=hx, hz=hz, name="t", meta={})
    if code.k < kmin or code.components() > 1:
        return []
    if n_probe:
        hx, hz, _ = tile.build(x_tile, B, *probe_size(B, n_probe))
        code = CSSCode(hx=hx, hz=hz, name="t", meta={})
        if code.k < kmin or fr.screen_distance(code, w, frontier, trials) is None:
            return []
    out = []
    for l, m in bulk_sizes(B, n_max, aspect=aspect):
        hx, hz, _ = tile.build(x_tile, B, l, m)
        n = hx.shape[1]
        k = n - gf2.rank(hx) - gf2.rank(hz)
        if k < kmin:
            continue
        if fr.bar(n, k, w, frontier) > n:                      # nothing this size could be new
            continue
        code = CSSCode(hx=hx, hz=hz, name="t", meta={})
        d = fr.screen_distance(code, w, frontier, trials)      # stops at the first logical below the bar
        if d is None:
            continue
        out.append({"tile": [list(e) for e in sorted(x_tile)], "B": B, "l": l, "m": m, "n": n, "k": k, "d_ub": d,
                    "w": w, "radius": round(tile_radius(x_tile, B), 3), "key": tile_key(x_tile, B),
                    "sig": code.signature()})
    return out


def cmd_screen(a) -> None:
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    frontier = fr.load("single", a.w)
    seen = {h for p in (LAB / "results" / "hunt").glob("*/screened.txt") for h in p.read_text().split()}
    tiles = [t for t in all_tiles(a.B, a.w, a.samples, a.seed) if key_hash(tile_key(t, a.B)) not in seen]
    print(f"{len(tiles)} tiles (B={a.B}, w={a.w}) to screen against {len(frontier)} frontier entries", flush=True)
    t_end, kept = time.time() + a.hours * 3600, 0
    with ProcessPoolExecutor(a.workers) as ex, open(out / "screen.jsonl", "a") as f, open(out / "screened.txt", "a") as fs:
        futs = {ex.submit(screen_one, t, a.B, a.w, frontier, a.kmin, a.n_max, a.trials, a.cap, a.n_probe, a.aspect):
                tile_key(t, a.B)
                for t in tiles}
        for fut in as_completed(futs):
            fs.write(key_hash(futs[fut]) + "\n")
            for r in fut.result():
                kept += 1
                f.write(json.dumps(r) + "\n")
                f.flush()
                print(f"keep [[{r['n']},{r['k']},<={r['d_ub']}]] w={r['w']} r={r['radius']} {r['l']}x{r['m']} {r['key']}",
                      flush=True)
            if time.time() > t_end:
                for x in futs:
                    x.cancel()
                break
    print(f"screen done: {kept} kept", flush=True)


def refine_one(r: dict, trials: int) -> dict | None:
    """Tighter d_ub before an exact proof: the screen's 100 trials overstate d at n > 300 (30 Sep: RIS said <= 19,
    DistQLDPC found 15). Uses any lighter logical an earlier proof run found; None once the frontier's bar is missed."""
    need = fr.bar(r["n"], r["k"], r["w"], fr.load("single", r["w"]))
    known = [s["ub"] for s in r.get("sides", {}).values() if s.get("ub")]
    d = min([r["d_ub"]] + known)
    if d >= need:
        hx, hz, _ = tile.build([tuple(e) for e in r["tile"]], r["B"], r["l"], r["m"])
        d = min(d, CSSCode(hx=hx, hz=hz, name="t", meta={}).distance_upper_bound(trials=trials, seed=7, stop_at=need - 1))
    return None if d < need else {**r, "d_ub": d, "refined": trials}


def cmd_prove(a) -> None:
    from qec_search.distqldpc import solve_many
    out = pathlib.Path(a.out)
    (out / "cand").mkdir(parents=True, exist_ok=True)
    last = {}
    for r in hb.read_jsonl(out / "proofs.jsonl"):
        last[(r["key"], r["l"], r["m"])] = r                  # a later line supersedes an earlier one
    done = {key for key, r in last.items() if not (a.retry_open and r["d_exact"] is None)}
    todo, seen_sig = [], set()
    for r in sorted(hb.read_jsonl(out / "screen.jsonl"), key=lambda r: -r["k"] * r["d_ub"] ** 2 / r["n"]):
        if (r["key"], r["l"], r["m"]) not in done and r["sig"] not in seen_sig:
            seen_sig.add(r["sig"])
            todo.append(r)
    todo = [{**r, **({"sides": last[(r["key"], r["l"], r["m"])]["sides"]} if (r["key"], r["l"], r["m"]) in last else {})}
            for r in todo[:3 * a.max]]
    with ProcessPoolExecutor(a.workers) as ex:
        todo = [r for r in ex.map(refine_one, todo, [a.refine_trials] * len(todo)) if r is not None]
    todo = sorted(todo, key=lambda r: -r["k"] * r["d_ub"] ** 2 / r["n"])[:a.max]
    print(f"{len(todo)} codes to prove after refining", flush=True)
    codes = {}
    for r in todo:
        hx, hz, q = tile.build([tuple(e) for e in r["tile"]], r["B"], r["l"], r["m"])
        stem = f"tile_{r['n']}_{r['k']}_{key_hash(r['key'])}_{r['l']}x{r['m']}"
        np.savez(out / "cand" / f"{stem}.npz", hx=hx, hz=hz, coords=tile.layout(q))
        codes[stem] = (hx.astype(np.int8), hz.astype(np.int8), r)
    sol = solve_many({s: (hx, hz) for s, (hx, hz, _) in codes.items()}, a.cube_secs, a.jobs)
    with open(out / "proofs.jsonl", "a") as f:
        for stem, (hx, hz, r) in codes.items():
            exact = [s["exact"] for s in sol[stem]["sides"].values()]
            d = None if None in exact else min(exact)
            frontier = fr.load("single", r["w"])
            p = {**r, "npz": f"{out.as_posix()}/cand/{stem}.npz", "sides": sol[stem]["sides"], "d_exact": d,
                 "board_twin": hb.board_twin(r["n"], r["sig"]),
                 "new_entry": d is not None and not fr.dominated((r["n"], r["k"], d, r["w"]), frontier)}
            f.write(json.dumps(p) + "\n")
            print(f"[[{p['n']},{p['k']},{d}]] (ub {p['d_ub']}) new_entry={p['new_entry']} twin={p['board_twin']} {stem}",
                  flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("screen")
    s.add_argument("--B", type=int, default=3)
    s.add_argument("--w", type=int, default=6, choices=[4, 5, 6, 7, 8])
    s.add_argument("--cap", type=float, default=4.0, help="interaction radius cap (4.0 = single layer)")
    s.add_argument("--samples", type=int, default=50_000, help="random tiles when B > 3")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--kmin", type=int, default=2)
    s.add_argument("--n-max", type=int, default=400)
    s.add_argument("--trials", type=int, default=100)
    s.add_argument("--hours", type=float, default=4.0)
    s.add_argument("--n-probe", type=int, default=200, help="screen each tile at this n first (0: sweep every tile)")
    s.add_argument("--aspect", type=float, default=2.0, help="largest m / l in the size sweep")
    p = sub.add_parser("prove")
    p.add_argument("--max", type=int, default=20)
    p.add_argument("--cube-secs", type=int, default=1800)
    p.add_argument("--jobs", type=int, default=os.cpu_count())
    p.add_argument("--retry-open", action="store_true", help="re-prove codes whose last proof did not finish")
    p.add_argument("--refine-trials", type=int, default=2000, help="d_ub trials on the leaders before proving")
    for x in (s, p):
        x.add_argument("--out", required=True)
        x.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args(argv)
    {"screen": cmd_screen, "prove": cmd_prove}[a.cmd](a)


if __name__ == "__main__":
    main()
