# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Exact distances for board codes that carry only a witnessed upper bound (no file in certs/).

Each side is posed for DistQLDPC (scripts/export_distqldpc.py, one CSS side, complete logical basis) and solved with
a time cap. A side that finishes gives its exact distance. If the smaller side is below the board's d, the board
overclaims; a lighter logical is then fetched with certify_sym's CryptoMiniSat encoding at weight <= that value and
checked (in ker H_opp, anticommutes with a logical), because a distance revision needs an explicit witness.

Results go to results/logs/board_distance_audit.jsonl, one line per code, and a rerun skips codes already there, so
nothing is solved twice. Run from lab/ in Git Bash (MSYS_NO_PATHCONV=1 for the WSL call is set here):

    python scripts/audit_board_distances.py --local-only --n-max 400 --secs 900
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time

LAB = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB))
from qec_search.certify_sym import is_logical, load_code, solve_cube  # noqa: E402
from qec_search.distqldpc import parse_log, solve_many  # noqa: E402,F401  (parse_log re-exported for tests)

OUT = LAB / "results" / "logs" / "board_distance_audit.jsonl"


def verdict(claimed: int, sides: dict) -> str:
    """LOWER (some side has a logical lighter than the board's d), holds (both sides proven >= d), open (unfinished),
    or ERROR (both exact and above d: impossible, since the board's d comes with a witness)."""
    ub = [s["exact"] if s["exact"] is not None else s["ub"] for s in sides.values()]
    lb = [s["exact"] if s["exact"] is not None else (s["lb"] or 0) for s in sides.values()]
    if any(u is not None and u < claimed for u in ub):
        return "LOWER"
    if min(lb) < claimed:
        return "open"
    return "holds" if min(lb) == claimed or any(u == claimed for u in ub) else "ERROR"


def targets(board: pathlib.Path, local_only: bool, n_max: int, d_max: int) -> list[pathlib.Path]:
    certs = {p.stem for p in (board / "certs").glob("*.json")}
    local = None
    if local_only:
        entries = json.loads((LAB / "results" / "hunt" / "board_entries.json").read_text())
        local = {e["slug"] for e in entries if e["locality_class"] in ("local-2d-single", "local-2d-bilayer")}
    out = []
    for p in sorted(board.glob("codes/*.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        if d.get("code_type", "CSS") != "CSS" or p.stem in certs:
            continue
        if local is not None and p.stem not in local:
            continue
        if d["n"] <= n_max and d["distance"]["d"] <= d_max:
            out.append(p)
    return sorted(out, key=lambda p: json.loads(p.read_text(encoding="utf-8"))["n"])


def witness_below(code: pathlib.Path, side: str, w: int, secs: float) -> list[int] | None:
    res = solve_cube(str(code), side, w + 1, [], secs)             # asks for a logical of weight <= w
    return res.get("witness") if res["status"] == "SAT" else None


def audit_batch(codes: list[pathlib.Path], secs: int, jobs: int) -> list[dict]:
    """Exact distances through qec_search.distqldpc (duality, orbital cubes), all codes in one solver batch."""
    loaded = {p: load_code(p) for p in codes}
    sol = solve_many({f"board_{p.stem}": (hx, hz) for p, (hx, hz, _) in loaded.items()}, secs, jobs)
    rows = {p: sol[f"board_{p.stem}"]["sides"] for p in codes}
    out = []
    for p, sides in rows.items():
        doc = json.loads(p.read_text(encoding="utf-8"))
        claimed = doc["distance"]["d"]
        r = {"slug": p.stem, "n": doc["n"], "k": doc["k"], "claimed_d": claimed, "sides": sides,
             "verdict": verdict(claimed, sides), "secs_cap": secs, "method": "distqldpc+orbital",
             "date": time.strftime("%Y-%m-%d")}
        if r["verdict"] == "LOWER":
            side, s = min(sides.items(), key=lambda kv: kv[1]["exact"] if kv[1]["exact"] is not None else kv[1]["ub"])
            w = s["exact"] if s["exact"] is not None else s["ub"]
            wit = witness_below(p, side, w, 1800)
            hx, hz, _ = loaded[p]
            h_same, h_opp = (hx, hz) if side == "X" else (hz, hx)
            r["witness"] = {"side": side, "weight": len(wit) if wit else None, "support": wit,
                            "checked": bool(wit) and is_logical(h_same, h_opp, wit)}
        out.append(r)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--board", type=pathlib.Path, default=LAB.parent / ".worktrees" / "frontier",
                    help="an upstream checkout (codes/, certs/)")
    ap.add_argument("--local-only", action="store_true", help="only codes in a 2D-local cell")
    ap.add_argument("--n-max", type=int, default=400)
    ap.add_argument("--d-max", type=int, default=24)
    ap.add_argument("--secs", type=int, default=900, help="DistQLDPC cap per side")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--batch", type=int, default=8, help="codes per solver batch (2 sides each)")
    ap.add_argument("--hours", type=float, default=12.0)
    ap.add_argument("--out", type=pathlib.Path, default=OUT)
    ap.add_argument("--retry-open", action="store_true", help="re-solve codes whose last verdict was open")
    a = ap.parse_args()
    out_path = a.out
    last = {}
    if out_path.exists():
        for line in out_path.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                last[r["slug"]] = r["verdict"]                  # a later line supersedes an earlier one
    done = {s for s, v in last.items() if not (a.retry_open and v == "open")}
    todo = [p for p in targets(a.board, a.local_only, a.n_max, a.d_max) if p.stem not in done]
    print(f"{len(todo)} codes to audit ({len(done)} already in {out_path.name})", flush=True)
    t_end = time.time() + a.hours * 3600
    for i in range(0, len(todo), a.batch):
        if time.time() > t_end:
            break
        for r in audit_batch(todo[i:i + a.batch], a.secs, a.jobs):
            with open(out_path, "a") as f:
                f.write(json.dumps(r) + "\n")
            print(f"{r['slug']:>20} claimed {r['claimed_d']:>3}  {r['verdict']:6} "
                  + " ".join(f"{s}:{v['exact'] if v['exact'] is not None else (v['lb'], v['ub'])}"
                             for s, v in r["sides"].items()), flush=True)


if __name__ == "__main__":
    main()
