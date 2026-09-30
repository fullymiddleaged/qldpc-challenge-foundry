# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Exact distance certification for any CSS code, with orbital branching on its Tanner-graph automorphisms.

The question per side is upstream's (verify/sat_certify.py): is there a nontrivial logical of weight <= d - 1? The
encoding is also upstream's: native XOR clauses for the checks, one selector per logical generator, and a
sequential-counter weight bound. One deliberate difference: the logical generators come from upstream's
gf2.logical_basis, and a check confirms they span all k classes. Upstream's sat_certify._logicals used to take rows
rank(H_same):rank(R) of a row-reduced stack, which on most board codes spanned fewer than k classes, so its UNSAT did
not rule out a lighter logical. We reported it (issue #2273); upstream fixed it in #2343 (28 Sep 2026).

On top of that encoding:

  * Orbital branching (Ostrowski, Linderoth, Rossi, Smriglio, Math. Program. 126, 147 (2011)). Let the Tanner-graph
    automorphism group Aut have qubit orbits O_1..O_s with representatives r_i. Every lightest logical c has an
    image under Aut that contains some r_i. If c meets O_1 at all, it can be moved onto r_1. So the search splits
    into cubes {v[r_1]=1}, {v[r_2]=1, O_1 all 0}, {v[r_3]=1, O_1 and O_2 all 0}, and so on. Within cube i, the
    stabiliser of r_i preserves every constraint, so its orbits on the remaining free qubits branch again
    (--depth). The automorphisms are found from H alone (qec_search.automorphisms) and each one is checked against
    the check supports, so nothing depends on how the code was constructed.
  * X/Z duality. A permutation that swaps the X and Z check sets gives d_X = d_Z, so only the X side is solved.
  * Adaptive cubes. A cube that exceeds --tlim is split on its lowest free qubit and requeued. Every finished cube is
    written to a status file straight away, and a rerun resumes from that file.

UNSAT on every leaf proves no logical of weight < d on that side. SAT on any leaf refutes the claim, and the witness
is checked (in the kernel, of weight <= d - 1, anticommuting with some logical) before it is reported.

    python -m qec_search.certify_sym ../codes/180-20-14.json               # d from the file's claim
    python -m qec_search.certify_sym results/candidates/x.npz --d 16
    python -m qec_search.certify_sym ../codes/180-20-14.json --status      # progress only
Status: results/logs/<stem>.sym_d<d>.status.json
"""
from __future__ import annotations

import argparse
import datetime as dt
import importlib.util
import json
import os
import pathlib
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np

from .automorphisms import find_duality, qubit_orbits

LAB = pathlib.Path(__file__).resolve().parents[1]
ROOT = LAB.parent
LOGS = LAB / "results" / "logs"
REPLACE_RETRIES = 10


def _upstream():
    """Upstream's verify/sat_certify.py (for its gf2 module), loaded from the enclosing repo."""
    sys.path.insert(0, str(ROOT / "verify"))
    spec = importlib.util.spec_from_file_location("_qldpc_sat_certify", ROOT / "verify" / "sat_certify.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_code(path: pathlib.Path) -> tuple[np.ndarray, np.ndarray, int | None]:
    """(hx, hz, claimed d or None) from a board JSON or a candidate .npz."""
    if path.suffix == ".npz":
        z = np.load(path)
        return z["hx"].astype(np.int8), z["hz"].astype(np.int8), None
    doc = json.loads(path.read_text())
    n = doc["n"]

    def mat(rows: list[list[int]]) -> np.ndarray:
        h = np.zeros((len(rows), n), dtype=np.int8)
        for i, r in enumerate(rows):
            h[i, r] = 1
        return h

    return mat(doc["checks"]["X"]), mat(doc["checks"]["Z"]), doc.get("distance", {}).get("d")


def orbital_cubes(hx: np.ndarray, hz: np.ndarray, depth: int) -> list[list[int]]:
    """Disjoint cubes (1-based DIMACS unit literals) that together cover, up to automorphism, every nonzero v.

    depth 0 is a single empty cube (no symmetry used). Each further level branches on the orbits of the pointwise
    stabiliser of the representatives already fixed, and stops early once that stabiliser is trivial."""
    n = hx.shape[1]

    def expand(fixed: tuple[int, ...], zero: frozenset[int], level: int) -> list[list[int]]:
        base = [q + 1 for q in fixed] + [-(q + 1) for q in sorted(zero)]
        if level == depth:
            return [base]
        free = set(range(n)) - zero - set(fixed)
        orbits, gens = qubit_orbits(hx, hz, fixed=fixed, among=free)
        if not gens:                                    # trivial group: branching would only be a plain case split
            return [base]
        out, done = [], set()
        for orb in orbits:
            out += expand(fixed + (orb[0],), zero | done, level + 1)
            done |= set(orb)
        return out

    return expand((), frozenset(), 0)


def children(cube: list[int], n: int) -> list[list[int]]:
    """Split a cube on the lowest qubit it leaves free."""
    fixed = {abs(x) for x in cube}
    v = next(q for q in range(1, n + 1) if q not in fixed)
    return [cube + [v], cube + [-v]]


def cube_key(side: str, cube: list[int]) -> str:
    return f"{side}|" + ",".join(map(str, cube))


def leaves(st: dict) -> list[tuple[str, list[int], int, int]]:
    """Every current leaf as (side, cube, root index, depth below root); SPLIT cubes give way to their halves."""
    out, stack = [], [(s, c, i, 0) for s in st["sides"] for i, c in enumerate(st["roots"])]
    while stack:
        side, cube, root, depth = stack.pop()
        if st["cubes"].get(cube_key(side, cube), {}).get("status") == "SPLIT":
            stack += [(side, ch, root, depth + 1) for ch in children(cube, st["n"])]
        else:
            out.append((side, cube, root, depth))
    return out


def summarize(st: dict) -> dict:
    counts = {"UNSAT": 0, "SAT": 0, "pending": 0}
    roots_done = {(s, i): True for s in st["sides"] for i in range(len(st["roots"]))}
    for side, cube, root, _ in leaves(st):
        status = st["cubes"].get(cube_key(side, cube), {}).get("status", "pending")
        counts[status] += 1
        if status != "UNSAT":
            roots_done[(side, root)] = False
    verdict = "REFUTED" if counts["SAT"] else ("PROVEN" if counts["pending"] == 0 else "in progress")
    return {"verdict": verdict, "roots_done": f"{sum(roots_done.values())}/{len(roots_done)}",
            "leaves_unsat": counts["UNSAT"], "leaves_pending": counts["pending"],
            "splits": sum(c["status"] == "SPLIT" for c in st["cubes"].values()),
            "cpu_hours": round(sum(c["secs"] for c in st["cubes"].values()) / 3600, 3)}


def pending(st: dict) -> list[tuple[str, list[int]]]:
    if any(c["status"] == "SAT" for c in st["cubes"].values()):
        return []
    return [(s, c) for s, c, _, _ in leaves(st) if st["cubes"].get(cube_key(s, c), {}).get("status") != "UNSAT"]


def pairing_set(h_same: np.ndarray, h_opp: np.ndarray) -> np.ndarray:
    """The other side's logical generators: ker(h_same) modulo rowspace(h_opp), one row per logical class.

    A v in ker(h_opp) is a nontrivial logical exactly when it anticommutes with one of these rows, but only if the
    rows span all k classes. That is checked here, and a shortfall is an error, not a weaker proof."""
    gf2 = _upstream().gf2
    t = np.asarray(gf2.logical_basis(h_same, h_opp), dtype=np.int8)
    n = h_same.shape[1]
    k = n - gf2.rank(h_same) - gf2.rank(h_opp)
    spanned = gf2.rank(np.vstack([h_opp, t])) - gf2.rank(h_opp) if len(t) else 0
    if not (len(t) == spanned == k):
        raise RuntimeError(f"logical basis spans {spanned} of k={k} classes with {len(t)} rows")
    return t


def is_logical(h_same: np.ndarray, h_opp: np.ndarray, support: list[int]) -> bool:
    """support is in ker(h_opp) and outside rowspace(h_same)."""
    v = np.zeros(h_opp.shape[1], dtype=np.int64)
    v[support] = 1
    if (h_opp.astype(np.int64) @ v % 2).any():
        return False
    return bool((pairing_set(h_same, h_opp).astype(np.int64) @ v % 2).any())


def solve_cube(path: str, side: str, d: int, cube: list[int], tlim: float) -> dict:
    """One cube of one side at W = d - 1."""
    from pycryptosat import Solver
    from pysat.card import CardEnc, EncType

    hx, hz, _ = load_code(pathlib.Path(path))
    h_same, h_opp = (hx, hz) if side == "X" else (hz, hx)
    n = hx.shape[1]
    t0 = time.time()
    s = Solver(threads=1)
    for row in h_opp:
        s.add_xor_clause([int(j) + 1 for j in np.nonzero(row)[0]], False)
    top, selectors = n, []
    for t in pairing_set(h_same, h_opp):
        top += 1
        selectors.append(top)
        s.add_xor_clause([int(j) + 1 for j in np.nonzero(t)[0]] + [top], False)
    if not selectors:
        return {"status": "UNSAT", "secs": 0.0}
    s.add_clause(selectors)
    card = CardEnc.atmost(lits=list(range(1, n + 1)), bound=d - 1, top_id=top, encoding=EncType.seqcounter)
    for cl in card.clauses:
        s.add_clause(cl)
    for lit in cube:
        s.add_clause([lit])
    sat, sol = s.solve(time_limit=tlim)
    out = {"status": "TIMEOUT" if sat is None else ("SAT" if sat else "UNSAT"), "secs": round(time.time() - t0, 2)}
    if sat:
        w = sorted(j - 1 for j in range(1, n + 1) if sol[j])
        if len(w) > d - 1 or not is_logical(h_same, h_opp, w):
            raise RuntimeError(f"solver returned an invalid witness on {path} {side} {cube}")
        out["witness"] = w
    return out


def status_path(code: pathlib.Path, d: int, tag: str = "") -> pathlib.Path:
    return LOGS / f"{code.stem}.sym{tag}_d{d}.status.json"


def new_status(code: pathlib.Path, d: int, depth: int, use_duality: bool) -> dict:
    hx, hz, _ = load_code(code)
    t0 = time.time()
    roots = orbital_cubes(hx, hz, depth)
    dual = use_duality and find_duality(hx, hz) is not None
    orbits, _ = qubit_orbits(hx, hz)
    return {"code": str(code), "n": int(hx.shape[1]), "claimed_d": d, "depth": depth,
            "orbit_sizes": [len(o) for o in orbits], "xz_duality": dual, "sides": ["X"] if dual else ["X", "Z"],
            "roots": roots, "symmetry_secs": round(time.time() - t0, 2),
            "method": "upstream sat_certify selector encoding + orbital branching on Tanner-graph automorphisms, "
                      "adaptive cubes (lab/qec_search/certify_sym.py)",
            "started": dt.datetime.now().isoformat(timespec="seconds"), "cubes": {}}


def save_status(path: pathlib.Path, st: dict) -> None:
    st["updated"] = dt.datetime.now().isoformat(timespec="seconds")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(st, indent=1))
    for attempt in range(REPLACE_RETRIES):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:                         # Windows: a watcher or indexer briefly holds the target open
            if attempt == REPLACE_RETRIES - 1:
                raise
            time.sleep(0.2 * (attempt + 1))


def _keep_awake() -> None:
    """Ask Windows not to sleep while this process runs (the display may still turn off); released on exit."""
    if os.name == "nt":
        import ctypes
        ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001
        ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED)


def _lower_priority() -> None:
    if os.name == "nt":
        import ctypes
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
    else:
        os.nice(10)


def certify(codes: list[pathlib.Path], d_arg: int | None, depth: int, use_duality: bool, tlim: float, workers: int,
            tag: str = "", quiet: bool = False, max_cpu_hours: float | None = None) -> dict[pathlib.Path, dict]:
    """Run (or resume) every code to a verdict; returns the final status per code.

    With max_cpu_hours, a code that has used that much solver time stops splitting: its timed-out cubes stay pending
    (verdict "in progress"), and a later run with a larger budget resumes them."""
    states, paths = {}, {}
    for code in codes:
        d = d_arg or load_code(code)[2]
        if not d:
            raise SystemExit(f"{code}: no claimed distance; pass --d")
        paths[code] = status_path(code, d, tag)
        states[code] = (json.loads(paths[code].read_text()) if paths[code].exists()
                        else new_status(code, d, depth, use_duality))
        save_status(paths[code], states[code])

    def submit(ex, code, side, cube):
        return ex.submit(solve_cube, str(code), side, states[code]["claimed_d"], cube, tlim), (code, side, cube)

    def over_budget(st: dict) -> bool:
        return max_cpu_hours is not None and sum(c["secs"] for c in st["cubes"].values()) > max_cpu_hours * 3600

    with ProcessPoolExecutor(max_workers=workers, initializer=_lower_priority) as ex:
        futs = dict(submit(ex, c, s, cube) for c in codes if not over_budget(states[c]) for s, cube in pending(states[c]))
        if not quiet:
            print(f"{len(futs)} cubes queued on {workers} workers", flush=True)
        while futs:
            done, _ = wait(futs, return_when=FIRST_COMPLETED)
            for f in done:
                if f not in futs:                       # dropped when its code was refuted earlier in this batch
                    continue
                code, side, cube = futs.pop(f)
                st, res = states[code], f.result()
                if res["status"] == "TIMEOUT":
                    if over_budget(st):
                        continue                        # leave this cube pending
                    res["status"] = "SPLIT"
                    futs.update(submit(ex, code, side, ch) for ch in children(cube, st["n"]))
                st["cubes"][cube_key(side, cube)] = res
                summ = summarize(st)
                st["summary"] = summ
                save_status(paths[code], st)
                if not quiet:
                    print(f"{code.stem} {side} root-cube len {len(cube)} {res['status']} {res['secs']}s "
                          f"{json.dumps(summ)}", flush=True)
                if summ["verdict"] == "REFUTED":
                    for g in [g for g, job in futs.items() if job[0] == code]:
                        g.cancel()
                        futs.pop(g)
    return states


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("codes", nargs="+", type=pathlib.Path, help="board JSON files or candidate .npz files")
    ap.add_argument("--d", type=int, help="claimed distance (default: the JSON's distance.d)")
    ap.add_argument("--depth", type=int, default=2, help="orbital-branching levels (0 = no symmetry)")
    ap.add_argument("--no-duality", action="store_true", help="solve both sides even if an X/Z duality exists")
    ap.add_argument("--tlim", type=float, default=600, help="seconds per cube before it is split")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1, help="solver processes (below-normal priority)")
    ap.add_argument("--max-cpu-hours", type=float, help="stop splitting a code once it has used this much solver time")
    ap.add_argument("--tag", default="", help="suffix for the status file, to keep benchmark arms apart")
    ap.add_argument("--status", action="store_true", help="print progress and exit")
    args = ap.parse_args()
    if args.status:
        for code in args.codes:
            d = args.d or load_code(code)[2]
            p = status_path(code, d, args.tag)
            print(code.stem, json.dumps(summarize(json.loads(p.read_text())) if p.exists() else "not started"))
        return
    _keep_awake()
    states = certify(args.codes, args.d, args.depth, not args.no_duality, args.tlim, args.workers, args.tag,
                     max_cpu_hours=args.max_cpu_hours)
    for code, st in states.items():
        print(code.stem, json.dumps(summarize(st)))


if __name__ == "__main__":
    main()
