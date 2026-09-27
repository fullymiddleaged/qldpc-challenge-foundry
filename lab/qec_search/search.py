# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Evolutionary search over bivariate bicycle codes.

    python -m qec_search.search --hours 2 --out runs/run1

Every evaluated code with k > 0 is appended to <out>/results.jsonl as soon as
it is scored, so a run can be stopped at any time (Ctrl-C) and the results so
far are kept. Afterwards run

    python -m qec_search.summarize runs/run1

and send the report it writes back to Claude.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import random
import signal
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .bbcode import BBGenome, build_bb
from .decoders import backend_name
from .evaluate import code_capacity_ler
from .known_codes import KNOWN, KNOWN_SIGS, best_known_kd2n

# --------------------------------------------------------------------- config


def parse_sizes(s: str):
    out = []
    for part in s.split(","):
        l, m = part.lower().split("x")
        out.append((int(l), int(m)))
    return out


DEFAULT_SIZES = "6x6,9x6,12x6,15x3,8x8,10x6,12x4,14x6,15x5,16x6,12x8,18x6,21x3,10x10,12x10,24x6"


# ------------------------------------------------------------------ genomes


def random_term(l, m, mixed: bool, rng: random.Random):
    if mixed:
        return (rng.randrange(l), rng.randrange(m))
    return (rng.randrange(1, l), 0) if rng.random() < 0.5 else (0, rng.randrange(1, m))


def random_genome(sizes, ta, tb, mixed, rng) -> BBGenome:
    l, m = rng.choice(sizes)
    while True:
        # IBM-style: A mostly powers of one variable, B of the other; the
        # random_term mix keeps the door open to other shapes.
        A = tuple(random_term(l, m, mixed, rng) for _ in range(ta))
        B = tuple(random_term(l, m, mixed, rng) for _ in range(tb))
        g = BBGenome(l, m, A, B)
        if g.is_valid():
            return g.normalized()


def mutate(g: BBGenome, sizes, mixed, rng) -> BBGenome:
    for _ in range(50):
        A, B = list(g.A), list(g.B)
        l, m = g.l, g.m
        r = rng.random()
        if r < 0.15:
            # move to another torus size, keep exponents (mod new size)
            l, m = rng.choice(sizes)
        elif r < 0.25:
            A, B = B, A
        else:
            tgt = A if rng.random() < 0.5 else B
            i = rng.randrange(len(tgt))
            a, b = tgt[i]
            step = rng.choice([-2, -1, 1, 2])
            if rng.random() < 0.3:
                tgt[i] = random_term(l, m, mixed, rng)
            elif mixed:
                tgt[i] = ((a + step) % l, b) if rng.random() < 0.5 else (a, (b + step) % m)
            elif a:  # pure x-power: stay a pure x-power
                tgt[i] = ((a + step) % l, 0)
            else:  # pure y-power (or the constant 1)
                tgt[i] = (0, (b + step) % m)
        ng = BBGenome(l, m, tuple(A), tuple(B))
        if ng.is_valid():
            return ng.normalized()
    return random_genome(sizes, len(g.A), len(g.B), mixed, rng)


def crossover(a: BBGenome, b: BBGenome, rng) -> BBGenome:
    ng = BBGenome(a.l, a.m, a.A, b.B) if rng.random() < 0.5 else BBGenome(a.l, a.m, b.A, a.B)
    return ng.normalized() if ng.is_valid() else a


# --------------------------------------------------------------- evaluation


def evaluate_genome(gj: dict, cfg: dict) -> dict:
    """Runs in a worker process. Returns a JSON-serialisable record."""
    t0 = time.time()
    g = BBGenome.from_json(gj)
    code = build_bb(g)
    rec = {"key": g.canonical_key(), "genome": g.to_json(), "pretty": g.pretty(), "n": code.n, "k": code.k,
           "weight": len(g.A) + len(g.B)}
    if code.k == 0:
        rec["seconds"] = round(time.time() - t0, 3)
        return rec
    ncomp = code.components()
    if ncomp > 1:  # several independent copies of a smaller code
        rec["split"] = ncomp
        rec["seconds"] = round(time.time() - t0, 3)
        return rec
    rec["sig"] = code.signature()
    d = code.distance_upper_bound(trials=cfg["dist_trials"], seed=cfg["seed"])
    rec["d_ub"] = d
    rec["kd2n"] = round(code.k * d * d / code.n, 4)
    rec["rate"] = round(code.k / code.n, 4)
    ref = best_known_kd2n(code.n)
    rec["ref_kd2n"] = round(ref, 4)
    # Only spend Monte Carlo time on codes that are at least in the running.
    if d >= cfg["min_d"] and rec["kd2n"] >= cfg["ler_gate"] * ref:
        code.meta["d_ub"] = d
        # refine the distance bound for promising codes
        d2 = code.distance_upper_bound(trials=max(cfg["dist_trials"] * 4, code.n), seed=cfg["seed"] + 1)
        if d2 < d:
            rec["d_ub"] = d = d2
            rec["kd2n"] = round(code.k * d * d / code.n, 4)
        rec["cc"] = [code_capacity_ler(code, p, shots=cfg["shots"], max_fails=cfg["max_fails"],
                                       seed=cfg["seed"], force_numpy=cfg["force_numpy"])
                     for p in cfg["ps"]]
    rec["known"] = KNOWN_SIGS.get(rec["sig"])
    rec["seconds"] = round(time.time() - t0, 3)
    return rec


def fitness(rec: dict) -> float:
    """Selection score within a torus-size niche.

    k*d^2/n is the standard figure of merit (it is what makes the gross code
    look good). Measured decoding performance adds a bonus, using the
    pessimistic end of the confidence interval so a lucky zero-failure run
    with few shots is not over-rewarded.
    """
    if not rec.get("k") or rec.get("split") or rec.get("dup_of"):
        return -1.0
    f = rec.get("kd2n", 0.0)
    cc = rec.get("cc")
    if cc:
        # use the harshest noise level: at mild noise big codes all show 0
        # failures and the bonus cannot tell them apart
        hi = cc[-1]["ler_hi"]
        k = max(rec["k"], 1)
        lq = 1 - (1 - hi) ** (1 / k) if hi < 1 else 1.0
        f += 0.5 * (-math.log10(max(lq, 1e-6)))
    return f


def size_of(rec_or_genome) -> str:
    g = rec_or_genome["genome"] if isinstance(rec_or_genome, dict) else rec_or_genome.to_json()
    return f"{g['l']}x{g['m']}"


# ---------------------------------------------------------------------- main


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/run1")
    ap.add_argument("--hours", type=float, default=1.0, help="wall-clock budget")
    ap.add_argument("--max-evals", type=int, default=0, help="stop after this many evaluations (0 = no limit)")
    ap.add_argument("--sizes", default=DEFAULT_SIZES, help="comma list of torus sizes lxm")
    ap.add_argument("--terms", default="3,3", help="terms in A,B (IBM codes: 3,3)")
    ap.add_argument("--mixed", action="store_true", help="allow mixed monomials x^i y^j")
    ap.add_argument("--niche", type=int, default=12, help="population kept per torus size")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--ps", default="0.04,0.06", help="code-capacity noise levels")
    ap.add_argument("--shots", type=int, default=3000)
    ap.add_argument("--max-fails", type=int, default=150)
    ap.add_argument("--dist-trials", type=int, default=60)
    ap.add_argument("--min-d", type=int, default=4)
    ap.add_argument("--ler-gate", type=float, default=0.6,
                    help="run Monte Carlo only if k*d^2/n >= gate * best known at that n")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--seed-known", action="store_true", help="also start from the IBM codes whose sizes are searched")
    ap.add_argument("--force-numpy", action="store_true", help="use the slow fallback decoder")
    args = ap.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    res_path = out / "results.jsonl"
    zk_path = out / "zero_k_keys.txt"
    rng = random.Random(args.seed)
    sizes = parse_sizes(args.sizes)
    ta, tb = (int(x) for x in args.terms.split(","))
    cfg = {"dist_trials": args.dist_trials, "seed": args.seed, "min_d": args.min_d,
           "ler_gate": args.ler_gate, "ps": [float(x) for x in args.ps.split(",")],
           "shots": args.shots, "max_fails": args.max_fails, "force_numpy": args.force_numpy}

    seen: dict[str, dict] = {}
    if res_path.exists():  # resume
        for line in res_path.open():
            try:
                r = json.loads(line)
                seen[r["key"]] = r
            except Exception:
                pass
    if zk_path.exists():
        for line in zk_path.open():
            seen.setdefault(line.strip(), {"k": 0})
    prev_meta = {}
    if (out / "meta.json").exists():
        try:
            prev_meta = json.load((out / "meta.json").open())
        except Exception:
            prev_meta = {}
    meta = {"started": time.strftime("%Y-%m-%d %H:%M:%S"), "args": vars(args), "decoder": backend_name(args.force_numpy),
            "python": sys.version.split()[0], "platform": platform.platform(), "resumed_with": len(seen)}
    try:
        import stim
        meta["stim"] = stim.__version__
    except Exception:
        meta["stim"] = None
    try:
        import ldpc
        meta["ldpc"] = getattr(ldpc, "__version__", "installed")
    except Exception:
        meta["ldpc"] = None
    with (out / "meta.json").open("w") as f:
        json.dump(meta, f, indent=2)

    print(f"decoder backend: {meta['decoder']}  | workers: {args.workers} | out: {out}")
    if meta["decoder"].startswith("numpy"):
        print("  (ldpc not found - using slow fallback decoder; `pip install ldpc` for real runs)")

    # one sub-population ("niche") per torus size keeps the search from
    # collapsing onto whichever size happens to score highest
    size_keys = [f"{l}x{m}" for l, m in sizes]
    pop: dict[str, list] = {s: [] for s in size_keys}
    for r in seen.values():
        if r.get("k") and size_of(r) in pop and fitness(r) >= 0:
            pop[size_of(r)].append(r)
    for s in pop:
        pop[s] = sorted(pop[s], key=fitness, reverse=True)[: args.niche]
    best_by_size = {s: (fitness(v[0]) if v else -1) for s, v in pop.items()}

    queue: list[BBGenome] = []
    if args.seed_known:
        queue += [g.normalized() for _, g, _ in KNOWN if (g.l, g.m) in sizes]
    while len(queue) < len(sizes) * 4:
        queue.append(random_genome(sizes, ta, tb, args.mixed, rng))

    deadline = time.time() + args.hours * 3600
    evals = zero_k = dup = sig_dups = 0
    sig_seen: dict[str, str] = {}
    for r in seen.values():
        if r.get("sig") and not r.get("dup_of"):
            sig_seen.setdefault(r["sig"], r["key"])
    tried_by_size: dict[str, int] = dict(prev_meta.get("tried_by_size", {}))
    stop = False

    def handle_sigint(*_):
        nonlocal stop
        stop = True
        print("\nstopping after in-flight evaluations...")

    signal.signal(signal.SIGINT, handle_sigint)
    last_print = time.time()

    def next_child() -> BBGenome:
        filled = [s for s, v in pop.items() if len(v) >= 2]
        if filled and rng.random() < 0.85:
            s = rng.choice(filled)
            niche = pop[s]
            parents = [max(rng.sample(niche, min(3, len(niche))), key=fitness) for _ in range(2)]
            pa, pb = (BBGenome.from_json(p["genome"]) for p in parents)
            child = crossover(pa, pb, rng) if rng.random() < 0.3 else pa
            return mutate(child, sizes, args.mixed, rng)
        return random_genome(sizes, ta, tb, args.mixed, rng)

    with ProcessPoolExecutor(max_workers=args.workers) as ex, res_path.open("a") as fout, zk_path.open("a") as fzk:
        inflight = {}
        tries = 0
        while not stop and time.time() < deadline and (not args.max_evals or evals < args.max_evals):
            while len(inflight) < args.workers * 2 and tries < 5000:
                g = queue.pop() if queue else next_child()
                key = g.canonical_key()
                if key in seen or key in inflight.values():
                    dup += 1
                    tries += 1
                    continue
                tries = 0
                fut = ex.submit(evaluate_genome, g.to_json(), cfg)
                inflight[fut] = key
            if not inflight:
                print("search space around current population exhausted; stopping")
                break
            done_fut = next(as_completed(list(inflight)))
            key = inflight.pop(done_fut)
            rec = done_fut.result()
            rec["key"] = key
            evals += 1
            seen[key] = rec
            s = size_of(rec)
            tried_by_size[s] = tried_by_size.get(s, 0) + 1
            if not rec["k"]:
                zero_k += 1
                fzk.write(key + "\n")
                continue
            if rec.get("sig"):
                first = sig_seen.setdefault(rec["sig"], key)
                if first != key:
                    rec["dup_of"] = first
                    sig_dups += 1
            fout.write(json.dumps(rec) + "\n")
            fout.flush()
            if rec.get("split") or rec.get("dup_of"):
                continue
            if s in pop:
                pop[s].append(rec)
                pop[s].sort(key=fitness, reverse=True)
                del pop[s][args.niche:]
                f = fitness(rec)
                if f > best_by_size.get(s, -1):
                    best_by_size[s] = f
                    ler = f" LER@{cfg['ps'][0]}={rec['cc'][0]['ler']:.2e}" if rec.get("cc") else ""
                    print(f"[{evals}] best on {s}: [[{rec['n']},{rec['k']},{rec['d_ub']}]] k*d^2/n={rec['kd2n']}{ler}  {rec['pretty']}"
                          + (f"  (= {rec['known']})" if rec.get("known") else ""))
            if time.time() - last_print > 60:
                last_print = time.time()
                left = (deadline - time.time()) / 60
                print(f"[{evals}] evaluated ({zero_k} had k=0, {dup + sig_dups} duplicates skipped), {left:.0f} min left")
        for fut in inflight:
            fut.cancel()

    meta["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    meta["evals_total"] = prev_meta.get("evals_total", 0) + evals
    meta["zero_k_total"] = prev_meta.get("zero_k_total", 0) + zero_k
    meta["evals_this_run"] = evals
    meta["zero_k_this_run"] = zero_k
    meta["duplicates_skipped"] = dup
    meta["relabelled_duplicates"] = sig_dups
    meta["tried_by_size"] = tried_by_size
    with (out / "meta.json").open("w") as f:
        json.dump(meta, f, indent=2)
    print(f"done: {evals} evaluations ({zero_k} with k=0). Now run:  python -m qec_search.summarize {out}")


if __name__ == "__main__":
    main()
