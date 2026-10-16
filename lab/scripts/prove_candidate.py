# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Exact distance for manifest candidates by DistQLDPC (qec_search.distqldpc.solve_many: each CSS side posed alone with a
complete logical basis, orbital cubes), for codes given only as an .npz (grafted tiles, anything without a genome).
Writes results/proofs/distance_proof_<n>_<k>_<d>.json per candidate, proven or not, in the shape submit.proof_refusal
reads; a side that did not finish within --secs leaves the result "open" with its proven bounds.

    python scripts/prove_candidate.py tile_360_8_17_a637_graft tile_404_8_18_a637_graft --secs 21600
--from-logs rebuilds the proof files from the solver logs without solving (the WSL solvers outlive a killed Windows
parent, so a run whose caller died can still be read back once its logs end with a "wall" line).
"""
from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import sys

import numpy as np

LAB = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB))
METHOD = ("DistQLDPC (MaxCDCL MaxSAT, arXiv:2606.12445) via qec_search.distqldpc: each CSS side posed alone with a "
          "complete logical basis; optimum per side (both sides)")


def params(info: dict) -> tuple[int, int, int]:
    n, k, d = (int(x) for x in info["params"].strip("[]").split(","))
    return n, k, d


def proof_record(name: str, info: dict, sides: dict, secs: int, today: str) -> dict:
    """The proof file for one candidate from solve_many's per-side results. Exact only when every side finished and
    the minimum over sides equals the manifest's claimed d; otherwise the record says what was proven."""
    n, k, d = params(info)
    exact = [s.get("exact") for s in sides.values()]
    if len(sides) == 2 and None not in exact:
        got = min(exact)
        result = "minimum distance proven exact" if got == d else f"proven d = {got}, manifest claims {d}"
    else:
        lb = min((s["exact"] if s.get("exact") is not None else (s.get("lb") or 0)) for s in sides.values()) if sides else 0
        result = f"open: proven d >= {lb} within {secs} s per cube"
    return {"code": f"[[{n},{k},{d}]]", "result": result, "d": d, "method": METHOD, "sides": sides, "date": today,
            "genome": info.get("genome"), "candidate": name,
            "reproduce": f"np.load('results/candidates/{name}.npz') (hx, hz) then qec_search.distqldpc.solve_many"}


def sides_from_logs(name: str, hx: np.ndarray, hz: np.ndarray) -> dict:
    """Per-side results from the logs of a previous solve_many run on this code."""
    from qec_search import distqldpc as dq
    per_side: dict[str, list[dict]] = {}
    for j in dq.plan(name, hx, hz):
        log = dq.DQ / "logs" / f"{j.stem}.log"
        per_side.setdefault(j.side, []).append(dq.parse_log(log.read_text() if log.exists() else ""))
    return {side: dq.combine(rs) for side, rs in per_side.items()}


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("names", nargs="+")
    p.add_argument("--secs", type=int, default=21600, help="time limit per cube")
    p.add_argument("--manifest", default=str(LAB / "results/candidates/manifest.json"))
    p.add_argument("--from-logs", action="store_true", help="read the solver logs of an earlier run; do not solve")
    a = p.parse_args(argv)
    from qec_search.distqldpc import solve_many

    man = json.loads(pathlib.Path(a.manifest).read_text(encoding="utf-8"))
    codes = {}
    for name in a.names:
        with np.load(LAB / "results/candidates" / f"{name}.npz") as z:
            codes[name] = (z["hx"].astype(np.int8), z["hz"].astype(np.int8))
    if a.from_logs:
        sol = {name: {"sides": sides_from_logs(name, *codes[name])} for name in a.names}
    else:
        sol = solve_many(codes, a.secs)
    today = datetime.date.today().isoformat()
    for name in a.names:
        rec = proof_record(name, man[name], sol[name]["sides"], a.secs, today)
        n, k, d = params(man[name])
        path = LAB / "results/proofs" / f"distance_proof_{n}_{k}_{d}.json"
        path.write_text(json.dumps(rec, indent=1) + "\n", encoding="utf-8")
        print(f"{name}: {rec['result']}  sides={rec['sides']}  -> {path.name}", flush=True)


if __name__ == "__main__":
    main()
