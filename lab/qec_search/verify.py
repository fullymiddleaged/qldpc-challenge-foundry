"""Slow loop: circuit-level evaluation of the best candidates (needs stim + ldpc).

    python -m qec_search.verify runs/run1 --top 5 --ps 0.003,0.002

Takes the top distinct [[n,k,d]] candidates from the run, plus the IBM
reference code of the same size as a control, builds a noisy syndrome
circuit for each, and appends results to <run>/circuit.jsonl. Re-run
summarize afterwards; the report gains a circuit-level table.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from .bbcode import BBGenome, build_bb
from .evaluate import circuit_level_ler
from .known_codes import KNOWN


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--ps", default="0.003")
    ap.add_argument("--rounds", type=int, default=0, help="0 = use the distance bound")
    ap.add_argument("--shots", type=int, default=5000)
    ap.add_argument("--max-fails", type=int, default=100)
    ap.add_argument("--no-controls", action="store_true", help="skip the IBM reference codes")
    args = ap.parse_args(argv)

    run = Path(args.run)
    tops = json.load((run / "top_candidates.json").open())
    picked, seen_params = [], set()
    for t in tops:
        if t["params"] in seen_params:
            continue
        seen_params.add(t["params"])
        picked.append(t)
        if len(picked) >= args.top:
            break
    todo = [(t["key"], BBGenome.from_json(t["genome"]), t["params"], t.get("known")) for t in picked]
    if not args.no_controls:
        ns = {int(t["params"][2:].split(",")[0]) for t in picked}
        for name, g, (n, k, d) in KNOWN:
            key = g.canonical_key()
            if n in ns and all(key != x[0] for x in todo):
                todo.append((key, g, f"{name} (control)", name))
    with (run / "circuit.jsonl").open("a") as f:
        for key, g, label, _known in todo:
            code = build_bb(g)
            d = code.distance_upper_bound(trials=200)
            code.meta["d_ub"] = d
            for p in (float(x) for x in args.ps.split(",")):
                t0 = time.time()
                res = circuit_level_ler(code, p, rounds=args.rounds or None, shots=args.shots, max_fails=args.max_fails)
                print(f"{label:40s} p={p} rounds={res['rounds']} LER={res['ler']:.2e} "
                      f"({res['fails']}/{res['shots']}) per-round={res['ler_per_round']:.2e}  {time.time() - t0:.0f}s")
                f.write(json.dumps({"key": key, "label": label, "genome": g.to_json(), "n": code.n, "k": code.k,
                                    "d_ub": d, "result": res}) + "\n")
                f.flush()
    # controls are not in results.jsonl; add minimal records so summarize can show them
    res_keys = {json.loads(l)["key"] for l in (run / "results.jsonl").open()} if (run / "results.jsonl").exists() else set()
    with (run / "results.jsonl").open("a") as f:
        for key, g, label, known in todo:
            if key not in res_keys:
                code = build_bb(g)
                d = code.distance_upper_bound(trials=200)
                f.write(json.dumps({"key": key, "genome": g.to_json(), "pretty": g.pretty(), "n": code.n, "k": code.k,
                                    "d_ub": d, "kd2n": round(code.k * d * d / code.n, 4),
                                    "known": known}) + "\n")
    print(f"done. Now re-run:  python -m qec_search.summarize {run}")


if __name__ == "__main__":
    main()
