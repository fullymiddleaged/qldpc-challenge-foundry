# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Find a 2D bilayer layout for a code so it qualifies for the qLDPC
Challenge's "2D-local bilayer" cells.

Rules we must satisfy (from the challenge's TRACKS.md / SCHEMA.md):
  * one [x, y] coordinate per qubit; `layers` = 2
  * at most 2 qubits share a site; distinct sites are >= 1 apart
    (we use an integer grid, so this holds automatically)
  * interaction radius = the largest *check diameter*, i.e. the largest
    distance between any two qubits in the same check. Bilayer cap: 7.0.

Method: start from a "folded torus" layout (each torus cell becomes one
site holding its left- and right-block qubit; each cyclic direction is
folded so wrap-around neighbours stay close), then simulated annealing
that swaps qubits between sites to push every check diameter under the
cap. Pure numpy; a few minutes per code on one CPU.

    python -m qec_search.layout --genome '{"l":18,"m":6,"A":[[0,5],[1,0],[12,0]],"B":[[0,1],[2,0],[3,0]]}'
"""
from __future__ import annotations

import argparse
import json
import math
import time

import numpy as np

from .bbcode import BBGenome, CSSCode, build_bb


def check_supports(code: CSSCode) -> np.ndarray:
    rows = [np.nonzero(r)[0] for r in np.vstack([code.hx, code.hz])]
    w = max(len(r) for r in rows)
    out = np.full((len(rows), w), -1, dtype=np.int64)
    for i, r in enumerate(rows):
        out[i, : len(r)] = r
        out[i, len(r):] = r[0]  # pad with a repeat; does not change the diameter
    return out


def diameters(coords: np.ndarray, supp: np.ndarray) -> np.ndarray:
    p = coords[supp]  # (c, w, 2)
    d = p[:, :, None, :] - p[:, None, :, :]
    return np.sqrt((d * d).sum(-1)).max(axis=(1, 2))


def fold(i: int, L: int) -> int:
    """Embed the cycle 0..L-1 in a line so cycle neighbours are <= 2 apart."""
    h = (L + 1) // 2
    return 2 * i if i < h else 2 * (L - i) - 1


def folded_torus(l: int, m: int, fold_x=True, fold_y=True) -> np.ndarray:
    n = 2 * l * m
    coords = np.zeros((n, 2))
    for i in range(l):
        for j in range(m):
            c = (fold(i, l) if fold_x else i, fold(j, m) if fold_y else j)
            coords[i * m + j] = c
            coords[l * m + i * m + j] = c
    return coords


def validate(coords: np.ndarray, supp: np.ndarray, layers: int = 2) -> dict:
    sites, counts = np.unique(np.round(coords, 9), axis=0, return_counts=True)
    dmin = math.inf
    if len(sites) > 1:
        diff = sites[:, None, :] - sites[None, :, :]
        dist = np.sqrt((diff ** 2).sum(-1)) + np.eye(len(sites)) * 1e9
        dmin = float(dist.min())
    diam = diameters(coords, supp)
    return {"max_per_site": int(counts.max()), "min_site_distance": dmin,
            "radius": float(diam.max()), "ok_sites": counts.max() <= layers and dmin >= 1 - 1e-9}


def anneal(code: CSSCode, init: np.ndarray, target: float = 7.0, width: int | None = None,
           height: int | None = None, steps: int = 400_000, seed: int = 0, verbose: bool = False):
    """Simulated annealing over a site grid with capacity 2 per site."""
    rng = np.random.default_rng(seed)
    n = code.n
    supp = check_supports(code)
    q2c = [[] for _ in range(n)]
    for c, row in enumerate(supp):
        for q in set(row.tolist()):
            q2c[q].append(c)
    q2c = [np.array(x) for x in q2c]

    xs, ys = init[:, 0].astype(int), init[:, 1].astype(int)
    W = width or int(xs.max()) + 1
    H = height or int(ys.max()) + 1
    # site occupancy: slot lists; qubits are placed into sites (cap 2)
    site_of = (xs * H + ys).astype(int)
    cap = np.zeros(W * H, dtype=int)
    for s in site_of:
        cap[s] += 1
    assert cap.max() <= 2
    site_xy = np.array([(s // H, s % H) for s in range(W * H)], dtype=float)
    coords = site_xy[site_of].copy()
    diam = diameters(coords, supp)

    def pen(d):
        over = np.maximum(d - target, 0)
        # any check over the cap costs at least 1, so the optimiser cannot
        # settle for "slightly over"
        return (over ** 2).sum() * 10 + (over > 0).sum() * 1.0 + (d ** 2).sum() * 1e-3

    cost = pen(diam)
    best = (diam.max(), coords.copy())
    T0, T1 = 2.0, 0.01
    t_start = time.time()
    for step in range(steps):
        T = T0 * (T1 / T0) ** (step / steps)
        a = int(rng.integers(n))
        if rng.random() < 0.5:
            b = int(rng.integers(n))  # swap two qubits
            if site_of[a] == site_of[b]:
                continue
            new_sa, new_sb = site_of[b], site_of[a]
            moved = [a, b]
        else:
            # move a to a nearby site with free capacity
            x, y = site_xy[site_of[a]]
            nx = int(np.clip(x + rng.integers(-3, 4), 0, W - 1))
            ny = int(np.clip(y + rng.integers(-3, 4), 0, H - 1))
            s = nx * H + ny
            if s == site_of[a] or cap[s] >= 2:
                continue
            new_sa, b, moved = s, None, [a]
        aff = np.unique(np.concatenate([q2c[q] for q in moved]))
        old = diam[aff]
        old_c = {q: coords[q].copy() for q in moved}
        coords[a] = site_xy[new_sa]
        if b is not None:
            coords[b] = site_xy[new_sb]
        new = diameters(coords, supp[aff])
        delta = pen(new) - pen(old)
        if delta <= 0 or rng.random() < math.exp(-delta / T):
            diam[aff] = new
            cost += delta
            if b is None:
                cap[site_of[a]] -= 1
                cap[new_sa] += 1
                site_of[a] = new_sa
            else:
                site_of[a], site_of[b] = new_sa, new_sb
            if diam.max() < best[0]:
                best = (diam.max(), coords.copy())
        else:
            for q in moved:
                coords[q] = old_c[q]
        if verbose and step % 50_000 == 0:
            print(f"  step {step:>7}  T={T:.3f}  max diameter {diam.max():.2f}  best {best[0]:.2f}  "
                  f"({time.time() - t_start:.0f}s)", flush=True)
        if best[0] <= target - 1e-9 and step > steps // 10:
            break
    return best[1], float(best[0])


def find_layout(code: CSSCode, genome: BBGenome, target: float = 7.0, steps: int = 400_000,
                restarts: int = 3, pad: int = 2, verbose: bool = False):
    supp = check_supports(code)
    best = (math.inf, None)
    starts = []
    for fx in (True, False):
        for fy in (True, False):
            starts.append((fx, fy))
    # score the structured starts, anneal from the best few
    scored = []
    for fx, fy in starts:
        c = folded_torus(genome.l, genome.m, fx, fy)
        scored.append((diameters(c, supp).max(), fx, fy, c))
    scored.sort(key=lambda t: t[0])
    for r in range(restarts):
        d0, fx, fy, c = scored[r % len(scored)]
        c = c + pad  # leave empty margin so qubits can spread out
        W = int(c[:, 0].max()) + 1 + pad
        H = int(c[:, 1].max()) + 1 + pad
        if verbose:
            print(f"restart {r}: start fold_x={fx} fold_y={fy}, grid {W}x{H}, start max diameter {d0:.2f}")
        coords, dmax = anneal(code, c, target=target, width=W, height=H, steps=steps, seed=r, verbose=verbose)
        if dmax < best[0]:
            best = (dmax, coords)
        if best[0] <= target:
            break
    return best[1], best[0]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True, help="genome JSON as in results.jsonl")
    ap.add_argument("--target", type=float, default=7.0)
    ap.add_argument("--steps", type=int, default=400_000)
    ap.add_argument("--restarts", type=int, default=3)
    ap.add_argument("--out", default=None, help="write hx, hz, coords to this .npz")
    args = ap.parse_args(argv)
    g = BBGenome.from_json(args.genome)
    code = build_bb(g)
    coords, dmax = find_layout(code, g, target=args.target, steps=args.steps, restarts=args.restarts, verbose=True)
    v = validate(coords, check_supports(code))
    print(f"result: max check diameter {v['radius']:.3f} (cap {args.target}), max qubits/site {v['max_per_site']}, "
          f"min site distance {v['min_site_distance']:.2f} -> {'QUALIFIES' if v['radius'] <= args.target and v['ok_sites'] else 'does not qualify'}")
    if args.out:
        np.savez(args.out, hx=code.hx.astype(np.uint8), hz=code.hz.astype(np.uint8), coords=coords)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
