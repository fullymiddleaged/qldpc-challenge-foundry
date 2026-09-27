"""Resumable, adaptive exact-distance certification of BB candidates, in the encoding of upstream's
verify/sat_certify.py.

Per side, the question is upstream's: is there a nontrivial logical of weight <= d - 1? It uses the same XOR
parity clauses, one selector per logical generator, and the sequential-counter weight bound. On top of that:
  * Torus translations. The script checks that the code is invariant under the x- and y-shifts of Z_l x Z_m. Those
    shifts act transitively on each block's cells, so a lightest logical can be moved to touch cell 0 of the left
    or right block. That gives the clause v[0] or v[lm], which is valid because translations preserve weight and
    triviality.
  * X/Z duality. If the relabelling (L,i,j) <-> (R,-i,-j) maps rowspace(H_Z) onto rowspace(H_X), which the script
    checks, it maps Z logicals onto X logicals of the same weight, so d_Z = d_X and only the X side is solved.
  * Adaptive cubes. Each side starts as disjoint cubes: v[0]=1, or v[0]=0 and v[lm]=1, times every assignment of
    --cube-bits further qubits. A cube that runs past --tlim is split in two on the next free qubit, and only its
    halves are requeued, so effort goes to the hard regions. Every finished cube is saved to the status file at
    once, so stopping loses at most one --tlim per worker, and a re-run resumes from the file.
UNSAT on every leaf cube proves d exactly; SAT on any cube refutes the claim and records the witness.

Run from the repo root so the solver deps are present:
    uv run --with pycryptosat --with python-sat python lab/scripts/certify_upstream.py \
        lab/results/candidates/bb_216_8_16_18x6.npz lab/results/candidates/bb_288_12_16_24x6.npz --d 16
    python lab/scripts/certify_upstream.py <candidates...> --d 16 --status     # progress only, no solving
Status: lab/results/logs/<candidate>.sat_d<d>.status.json
"""
from __future__ import annotations

import argparse
import datetime as dt
import importlib.util
import itertools
import json
import os
import pathlib
import re
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
LOGS = ROOT / "lab" / "results" / "logs"


def _sat_certify():
    sys.path.insert(0, str(ROOT / "verify"))
    spec = importlib.util.spec_from_file_location("_qldpc_sat_certify", ROOT / "verify" / "sat_certify.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def torus_shape(path: pathlib.Path) -> tuple[int, int]:
    """(l, m) from a candidate name such as bb_216_8_16_18x6."""
    match = re.search(r"_(\d+)x(\d+)$", path.stem)
    if not match:
        raise ValueError(f"no <l>x<m> torus in {path.name}")
    return int(match.group(1)), int(match.group(2))


def _cell(i: int, j: int, l: int, m: int) -> int:
    return (i % l) * m + (j % m)


def translation_perms(l: int, m: int) -> list[np.ndarray]:
    """Column permutations for the x- and y-shifts, acting on both blocks (qubit index = block*lm + i*m + j)."""
    perms = []
    for di, dj in ((1, 0), (0, 1)):
        cell = np.array([_cell(i + di, j + dj, l, m) for i in range(l) for j in range(m)])
        perms.append(np.concatenate([cell, cell + l * m]))
    return perms


def duality_perm(l: int, m: int) -> np.ndarray:
    """(L,i,j) <-> (R,-i,-j), an involution on the 2lm qubits."""
    cell = np.array([_cell(-i, -j, l, m) for i in range(l) for j in range(m)])
    return np.concatenate([cell + l * m, cell])


def _same_rowspace(a: np.ndarray, b: np.ndarray) -> bool:
    gf2 = _sat_certify().gf2
    ra = gf2.rank(a)
    return ra == gf2.rank(b) == gf2.rank(np.vstack([a, b]))


def has_translation_symmetry(hx: np.ndarray, hz: np.ndarray, l: int, m: int) -> bool:
    return all(_same_rowspace(h, h[:, p]) for p in translation_perms(l, m) for h in (hx, hz))


def has_xz_duality(hx: np.ndarray, hz: np.ndarray, l: int, m: int) -> bool:
    return _same_rowspace(hx, hz[:, duality_perm(l, m)])


def root_cubes(l: int, m: int, bits: int) -> list[list[int]]:
    """Disjoint unit-literal cubes (1-based DIMACS) covering every v with v[0] or v[lm] set."""
    heads = [[1], [-1, l * m + 1]]
    free = list(range(2, bits + 2))                     # qubits 1..bits of the left block
    return [head + [q if b else -q for q, b in zip(free, signs)]
            for head in heads for signs in itertools.product((1, 0), repeat=bits)]


def children(cube: list[int], n: int) -> list[list[int]]:
    """Split a cube on the lowest qubit it leaves free."""
    fixed = {abs(x) for x in cube}
    v = next(q for q in range(1, n + 1) if q not in fixed)
    return [cube + [v], cube + [-v]]


def cube_key(side: str, cube: list[int]) -> str:
    return f"{side}|" + ",".join(map(str, cube))


def leaves(st: dict, roots: list[list[int]], sides: list[str], n: int) -> list[tuple[str, list[int], int]]:
    """Every current leaf as (side, cube, depth below its root); SPLIT cubes are replaced by their halves."""
    out, stack = [], [(s, c, 0) for s in sides for c in roots]
    while stack:
        side, cube, depth = stack.pop()
        if st["cubes"].get(cube_key(side, cube), {}).get("status") == "SPLIT":
            stack += [(side, ch, depth + 1) for ch in children(cube, n)]
        else:
            out.append((side, cube, depth))
    return out


def summarize(st: dict, roots: list[list[int]], n: int) -> dict:
    """Counts of leaf cubes by status, and the share of the search space (weighted by cube size) proven empty."""
    sides = st["sides"]
    counts = {"UNSAT": 0, "SAT": 0, "pending": 0}
    proven = 0.0
    for side, cube, depth in leaves(st, roots, sides, n):
        status = st["cubes"].get(cube_key(side, cube), {}).get("status", "pending")
        counts[status] += 1
        if status == "UNSAT":
            proven += 0.5 ** depth / (len(roots) * len(sides))
    verdict = "REFUTED" if counts["SAT"] else ("PROVEN" if counts["pending"] == 0 else "in progress")
    return {"verdict": verdict, "proven_share": round(proven, 4), "leaves_unsat": counts["UNSAT"],
            "leaves_pending": counts["pending"], "splits": sum(c["status"] == "SPLIT" for c in st["cubes"].values()),
            "cpu_hours": round(sum(c["secs"] for c in st["cubes"].values()) / 3600, 2)}


def pending(st: dict, roots: list[list[int]], n: int) -> list[tuple[str, list[int]]]:
    """Leaves still to run; none once any cube is SAT, since the claim is then refuted."""
    if any(c["status"] == "SAT" for c in st["cubes"].values()):
        return []
    return [(s, c) for s, c, _ in leaves(st, roots, st["sides"], n)
            if st["cubes"].get(cube_key(s, c), {}).get("status") != "UNSAT"]


def _lower_priority() -> None:
    if os.name == "nt":
        import ctypes
        BELOW_NORMAL = 0x4000
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), BELOW_NORMAL)
    else:
        os.nice(10)


def solve_cube(path: str, side: str, d: int, cube: list[int], tlim: float) -> dict:
    """One cube of one side: upstream's selector encoding at W = d - 1, plus the translation clause and the cube."""
    from pycryptosat import Solver
    from pysat.card import CardEnc, EncType

    sc = _sat_certify()
    z = np.load(path)
    hx, hz = z["hx"].astype(np.int8), z["hz"].astype(np.int8)
    h_same, h_opp = (hx, hz) if side == "X" else (hz, hx)
    n = hx.shape[1]
    t0 = time.time()
    s = Solver(threads=1)
    for row in h_opp:
        s.add_xor_clause([int(j) + 1 for j in np.nonzero(row)[0]], False)
    top, selectors = n, []
    for t in sc._logicals(h_opp, h_same):
        top += 1
        selectors.append(top)
        s.add_xor_clause([int(j) + 1 for j in np.nonzero(t)[0]] + [top], False)
    s.add_clause(selectors)
    card = CardEnc.atmost(lits=list(range(1, n + 1)), bound=d - 1, top_id=top, encoding=EncType.seqcounter)
    for cl in card.clauses:
        s.add_clause(cl)
    s.add_clause([1, n // 2 + 1])                        # touches cell 0 of the left or right block
    for lit in cube:
        s.add_clause([lit])
    sat, sol = s.solve(time_limit=tlim)
    out = {"status": "TIMEOUT" if sat is None else ("SAT" if sat else "UNSAT"), "secs": round(time.time() - t0, 1)}
    if sat:
        out["witness"] = sorted(j - 1 for j in range(1, n + 1) if sol[j])
    return out


def status_path(cand: pathlib.Path, d: int) -> pathlib.Path:
    return LOGS / f"{cand.stem}.sat_d{d}.status.json"


def load_status(cand: pathlib.Path, d: int, bits: int) -> dict | None:
    p = status_path(cand, d)
    if not p.exists():
        return None
    st = json.loads(p.read_text())
    if st["cube_bits"] != bits:
        raise SystemExit(f"{p.name} was made with --cube-bits {st['cube_bits']}; rerun with that value")
    return st


def new_status(cand: pathlib.Path, d: int, bits: int) -> dict:
    """Check the symmetries on the matrices and record which sides need solving."""
    l, m = torus_shape(cand)
    z = np.load(cand)
    hx, hz = z["hx"].astype(np.int8), z["hz"].astype(np.int8)
    if not has_translation_symmetry(hx, hz, l, m):
        raise SystemExit(f"{cand.name}: not invariant under the Z_{l} x Z_{m} translations; the cell-0 clause is invalid")
    dual = has_xz_duality(hx, hz, l, m)
    return {"candidate": cand.name, "claimed_d": d, "cube_bits": bits, "sides": ["X"] if dual else ["X", "Z"],
            "xz_duality": dual, "translation_symmetry": True,
            "method": "upstream verify/sat_certify.py selector encoding + torus-translation clause, adaptive cubes",
            "cubes": {}}


def save_status(cand: pathlib.Path, d: int, st: dict) -> None:
    st["updated"] = dt.datetime.now().isoformat(timespec="seconds")
    p = status_path(cand, d)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(st, indent=1))
    os.replace(tmp, p)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates", nargs="+", type=pathlib.Path)
    ap.add_argument("--d", type=int, required=True, help="claimed distance on both sides")
    ap.add_argument("--tlim", type=float, default=600, help="seconds per cube before it is split")
    ap.add_argument("--cube-bits", type=int, default=2, help="root cubes per side = 2 * 2**bits")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help="solver processes (default: half the cores, at below-normal priority)")
    ap.add_argument("--status", action="store_true", help="print progress and exit")
    args = ap.parse_args()

    states, roots, sizes = {}, {}, {}
    for cand in args.candidates:
        l, m = torus_shape(cand)
        roots[cand], sizes[cand] = root_cubes(l, m, args.cube_bits), 2 * l * m
        st = load_status(cand, args.d, args.cube_bits)
        if args.status:
            print(cand.stem, json.dumps(summarize(st, roots[cand], sizes[cand]) if st else {"verdict": "not started"}))
        else:
            states[cand] = st or new_status(cand, args.d, args.cube_bits)
    if args.status:
        return

    def submit(ex, cand, side, cube):
        return ex.submit(solve_cube, str(cand), side, args.d, cube, args.tlim), (cand, side, cube)

    with ProcessPoolExecutor(max_workers=args.workers, initializer=_lower_priority) as ex:
        futs = dict(submit(ex, c, s, cube) for c in args.candidates for s, cube in pending(states[c], roots[c], sizes[c]))
        print(f"{len(futs)} cubes queued on {args.workers} workers", flush=True)
        while futs:
            done, _ = wait(futs, return_when=FIRST_COMPLETED)
            for f in done:
                cand, side, cube = futs.pop(f)
                st, res = states[cand], f.result()
                if res["status"] == "TIMEOUT":
                    res["status"] = "SPLIT"
                    futs.update(submit(ex, cand, side, ch) for ch in children(cube, sizes[cand]))
                st["cubes"][cube_key(side, cube)] = res
                save_status(cand, args.d, st)
                summ = summarize(st, roots[cand], sizes[cand])
                print(f"{cand.stem} {cube_key(side, cube)} {res['status']} {res['secs']}s  {json.dumps(summ)}",
                      flush=True)
                if summ["verdict"] == "REFUTED":            # drop the rest of this candidate's cubes
                    for g in [g for g, job in futs.items() if job[0] == cand]:
                        g.cancel()
                        futs.pop(g)
    for cand in args.candidates:
        print(cand.stem, json.dumps(summarize(states[cand], roots[cand], sizes[cand])))


if __name__ == "__main__":
    main()
