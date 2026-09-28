# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Hunt for bivariate bicycle codes that would be new nondominated entries in a 2D-local bilayer cell.

Three stages, each resumable from its output file in --out:
  screen  generate genomes (covers of given bases, or random), keep connected codes whose (n, k, d_ub, w) no board
          entry in the cell dominates. d_ub is an upper bound, so this only ever keeps too much, never too little.
  layout  anneal a bilayer layout (qec_search.layout) for each survivor, keep radius <= 7.0, write an .npz.
  prove   exact distance with DistQLDPC in WSL (scripts/run_distqldpc.sh), one job per orbital cube, then recheck
          dominance with the proven d.

Covers (arXiv:2511.13560): a term x^i y^j of a base code on Z_l x Z_m lifts to x^(i + t l) y^j on Z_(hl) x Z_m for
t in 0..h-1 (or likewise along y); every choice of t per term gives an h-fold cover of the base Tanner graph.

Each weight is judged in its own cell. The frontier for weight w is what `qldpc targets --cell "<cell>" --n 700
--top 400` prints for the cell "2D-local bilayer / weight <= w" (with the Unicode sign), saved as
results/hunt/frontier_bilayer_w<w>.txt; --frontier overrides it.

    python scripts/hunt_bilayer.py screen --out results/hunt/covers --covers bb_216_8_16_18x6,bb_288_12_16_24x6 --h 2
    python scripts/hunt_bilayer.py screen --out results/hunt/w8 --weight 8 --sizes 12x6,14x6,16x6 --hours 4
    python scripts/hunt_bilayer.py layout --out results/hunt/covers
    python scripts/hunt_bilayer.py prove --out results/hunt/covers --cube-secs 3600
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import pathlib
import random
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

LAB = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB))
from qec_search.bbcode import BBGenome, build_bb  # noqa: E402

CAP = 7.0
WSL_LAB = "/mnt/c/vscode/qldpc-challenge-foundry/lab"


# ------------------------------------------------------------------ frontier

def parse_frontier(text: str) -> list[tuple[int, int, int, int]]:
    """(n, k, d, w) for every '[[n,k,d]] w=W' entry in qldpc targets output."""
    return [tuple(map(int, m)) for m in re.findall(r"\[\[(\d+),(\d+),(\d+)\]\]\s+w=(\d+)", text)]


def load_frontier(a, w: int) -> list[tuple[int, int, int, int]]:
    p = pathlib.Path(a.frontier) if a.frontier else LAB / "results" / "hunt" / f"frontier_bilayer_w{w}.txt"
    return parse_frontier(p.read_text(encoding="utf-8"))


def dominated(c: tuple[int, int, int, int], frontier) -> bool:
    """Some entry is at least as good on n (lower), k, d (higher) and w (lower). Equal parameters count as dominated:
    a tie is no new record."""
    n, k, d, w = c
    return any(fn <= n and fk >= k and fd >= d and fw <= w for fn, fk, fd, fw in frontier)


# ------------------------------------------------------------------ genomes

def covers(g: BBGenome, h: int) -> list[BBGenome]:
    """All h-fold covers of g along x and along y, one per choice of lift for each term (deduplicated)."""
    out, seen = [], set()
    terms = list(g.A) + list(g.B)
    na = len(g.A)
    for axis in ("x", "y"):
        L, M = (h * g.l, g.m) if axis == "x" else (g.l, h * g.m)
        for lifts in itertools.product(range(h), repeat=len(terms)):
            new = [((i + t * g.l) % L, j) if axis == "x" else (i, (j + t * g.m) % M)
                   for (i, j), t in zip(terms, lifts)]
            c = BBGenome(L, M, tuple(new[:na]), tuple(new[na:]))
            if c.is_valid() and (key := c.canonical_key()) not in seen:
                seen.add(key)
                out.append(c.normalized())
    return out


def random_genomes(sizes, terms: int, count: int, seed: int) -> list[BBGenome]:
    from qec_search.search import random_genome
    rng = random.Random(seed)
    return [random_genome(sizes, terms, terms, False, rng) for _ in range(count)]


def stem_for(r: dict) -> str:
    g = BBGenome.from_json(r["genome"])
    return f"bb_{r['n']}_{r['k']}_{g.l}x{g.m}_{hashlib.sha1(r['key'].encode()).hexdigest()[:8]}"


# ------------------------------------------------------------------ stages

def screen_one(gj: dict, frontier, kmin: int, trials: int) -> dict | None:
    g = BBGenome.from_json(gj)
    code = build_bb(g)
    k = code.k
    if k < kmin or code.components() > 1:
        return None
    w = len(g.A) + len(g.B)
    d_ub = code.distance_upper_bound(trials=max(5, trials // 10))   # more trials only lower it, so a cheap
    if dominated((code.n, k, d_ub, w), frontier):                  # bound that is already beaten is final
        return None
    d_ub = min(d_ub, code.distance_upper_bound(trials=trials, seed=1))
    if dominated((code.n, k, d_ub, w), frontier):
        return None
    return {"genome": g.to_json(), "key": g.canonical_key(), "n": code.n, "k": k, "d_ub": d_ub, "w": w,
            "sig": code.signature()}


def read_jsonl(p: pathlib.Path) -> list[dict]:
    return [json.loads(line) for line in p.read_text().splitlines() if line.strip()] if p.exists() else []


def base_genomes(names: list[str]) -> list[BBGenome]:
    manifest = json.loads((LAB / "results" / "candidates" / "manifest.json").read_text())
    items = manifest if isinstance(manifest, list) else manifest.get("codes", manifest)
    items = [dict(name=k, **v) for k, v in items.items()] if isinstance(items, dict) else items
    found = {i["name"]: BBGenome.from_json(i["genome"]) for i in items if i["name"] in names}
    missing = set(names) - set(found)
    if missing:
        raise SystemExit(f"not in manifest: {sorted(missing)}")
    return list(found.values())


def cmd_screen(a) -> None:
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    frontier = load_frontier(a, 6 if a.covers else a.weight)
    done = {r["key"] for r in read_jsonl(out / "screen.jsonl")}
    if a.covers:
        gens = [c for b in base_genomes(a.covers.split(",")) for h in map(int, a.h.split(",")) for c in covers(b, h)]
    else:
        from qec_search.search import parse_sizes
        gens = random_genomes(parse_sizes(a.sizes), a.weight // 2, a.count, a.seed)
    uniq = {}
    for g in gens:
        if 2 * g.l * g.m <= a.n_max:
            uniq.setdefault(g.canonical_key(), g)
    gens = [g for key, g in uniq.items() if key not in done]
    print(f"{len(gens)} genomes to screen against {len(frontier)} frontier entries", flush=True)
    t_end = time.time() + a.hours * 3600
    kept = 0
    with ProcessPoolExecutor(a.workers) as ex, open(out / "screen.jsonl", "a") as f:
        futs = [ex.submit(screen_one, g.to_json(), frontier, a.kmin, a.trials) for g in gens]
        for fut in as_completed(futs):
            if (r := fut.result()) is not None:
                kept += 1
                f.write(json.dumps(r) + "\n")
                f.flush()
                print(f"keep [[{r['n']},{r['k']},<={r['d_ub']}]] w={r['w']} {r['key']}", flush=True)
            if time.time() > t_end:
                for x in futs:
                    x.cancel()
                break
    print(f"screen done: {kept} kept", flush=True)


def layout_one(r: dict, steps: int, restarts: int, npz_dir: str) -> dict:
    from qec_search.layout import check_supports, find_layout, validate
    g = BBGenome.from_json(r["genome"])
    code = build_bb(g)
    coords, _ = find_layout(code, g, target=CAP, steps=steps, restarts=restarts)
    v = validate(coords, check_supports(code))
    ok = v["radius"] <= CAP and v["ok_sites"]
    out = {**r, "radius": round(float(v["radius"]), 3), "layout_ok": bool(ok)}
    if ok:
        p = pathlib.Path(npz_dir) / f"{stem_for(r)}.npz"
        np.savez(p, hx=code.hx.astype(np.uint8), hz=code.hz.astype(np.uint8), coords=coords)
        out["npz"] = str(p.relative_to(LAB)).replace("\\", "/")
    return out


def cmd_layout(a) -> None:
    out = pathlib.Path(a.out)
    (out / "cand").mkdir(parents=True, exist_ok=True)
    done = {r["key"] for r in read_jsonl(out / "layout.jsonl")}
    todo, seen_sig = [], set()
    for r in sorted(read_jsonl(out / "screen.jsonl"), key=lambda r: -r["k"] * r["d_ub"] ** 2 / r["n"]):
        if r["key"] not in done and r["sig"] not in seen_sig:   # relabelled copies share a signature
            seen_sig.add(r["sig"])
            todo.append(r)
    todo = todo[:a.max]
    print(f"{len(todo)} layouts to anneal", flush=True)
    with ProcessPoolExecutor(a.workers) as ex, open(out / "layout.jsonl", "a") as f:
        for fut in as_completed([ex.submit(layout_one, r, a.steps, a.restarts, str((out / "cand").resolve()))
                                 for r in todo]):
            r = fut.result()
            f.write(json.dumps(r) + "\n")
            f.flush()
            print(f"[[{r['n']},{r['k']},<={r['d_ub']}]] radius {r['radius']} {'OK' if r['layout_ok'] else '-'}",
                  flush=True)


def prove_one(r: dict, frontier, cube_secs: int, jobs: int) -> dict:
    """Exact d over the orbital cubes of both sides (one side if an X/Z duality exists)."""
    from qec_search.automorphisms import find_duality
    from qec_search.certify_sym import orbital_cubes
    z = np.load(LAB / r["npz"])
    hx, hz = z["hx"].astype(np.int8), z["hz"].astype(np.int8)
    sides = ["X"] if find_duality(hx, hz) is not None else ["X", "Z"]
    stem = pathlib.Path(r["npz"]).stem
    args, tags = [], []
    for side in sides:
        for i, cube in enumerate(orbital_cubes(hx, hz, 2)):
            res = subprocess.run([sys.executable, str(LAB / "scripts" / "export_distqldpc.py"), str(LAB / r["npz"]),
                                  str(LAB / "results" / "distqldpc"), "--side", side,
                                  "--cube", ",".join(map(str, cube)), "--tag", f"c{i}"],
                                 capture_output=True, text=True, check=True).stdout.split()
            flag = next((t for t in res[1:] if t.startswith("-one-z=")), "")
            args.append(f"results/distqldpc/{stem}_{side}_c{i}" + (f":{flag}" if flag else ""))
            tags.append(f"{side}_c{i}")
    env = dict(os.environ, JOBS=str(jobs), WSLENV="JOBS")
    subprocess.run(["wsl.exe", "-d", "Ubuntu", "--", "bash", f"{WSL_LAB}/scripts/run_distqldpc.sh", str(cube_secs),
                    *args], env=env, capture_output=True)
    opt = {}
    for tag in tags:
        text = (LAB / "results" / "distqldpc" / "logs" / f"{stem}_{tag}.log").read_text()
        m = re.search(r"^o (\d+)", text, re.M)
        opt[tag] = int(m.group(1)) if m else None
    d = None if None in opt.values() else min(opt.values())
    res = {**r, "sides": sides, "cubes": opt, "d_exact": d}
    res["new_entry"] = d is not None and not dominated((r["n"], r["k"], d, r["w"]), frontier)
    return res


def cmd_prove(a) -> None:
    out = pathlib.Path(a.out)
    done = {r["key"] for r in read_jsonl(out / "proofs.jsonl")}
    todo = [r for r in read_jsonl(out / "layout.jsonl") if r["layout_ok"] and r["key"] not in done]
    todo.sort(key=lambda r: -r["k"] * r["d_ub"] ** 2 / r["n"])
    todo = todo[:a.max]
    print(f"{len(todo)} codes to prove", flush=True)
    with open(out / "proofs.jsonl", "a") as f:
        for r in todo:
            p = prove_one(r, load_frontier(a, r["w"]), a.cube_secs, a.jobs)
            f.write(json.dumps(p) + "\n")
            f.flush()
            print(f"[[{p['n']},{p['k']},{p['d_exact']}]] (ub {p['d_ub']}) new_entry={p['new_entry']} {p['npz']}",
                  flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("screen")
    s.add_argument("--covers", default="", help="manifest names of base codes to lift")
    s.add_argument("--h", default="2", help="cover degrees, comma-separated")
    s.add_argument("--sizes", default="12x6,14x6,16x6,18x6,20x6,22x6,24x6")
    s.add_argument("--weight", type=int, default=6, choices=[6, 8])
    s.add_argument("--count", type=int, default=200_000, help="random genomes to draw")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--kmin", type=int, default=4)
    s.add_argument("--n-max", type=int, default=700)
    s.add_argument("--trials", type=int, default=100, help="information-set trials for d_ub")
    s.add_argument("--hours", type=float, default=8.0)
    lay = sub.add_parser("layout")
    lay.add_argument("--steps", type=int, default=400_000)
    lay.add_argument("--restarts", type=int, default=2)
    p = sub.add_parser("prove")
    p.add_argument("--cube-secs", type=int, default=3600)
    p.add_argument("--jobs", type=int, default=os.cpu_count())
    for x in (lay, p):
        x.add_argument("--max", type=int, default=100, help="best codes by k d_ub^2 / n to take this run")
    for x in (s, lay, p):
        x.add_argument("--out", required=True)
        x.add_argument("--frontier", default="", help="override the per-weight frontier file")
        x.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args(argv)
    {"screen": cmd_screen, "layout": cmd_layout, "prove": cmd_prove}[a.cmd](a)


if __name__ == "__main__":
    main()
