# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Certify the minimum distance of a CSS code with an integer program.

The information-set search only ever finds logicals, so it gives an upper
bound. This module proves a LOWER bound: it asks a MILP solver (HiGHS,
bundled with scipy) for the minimum-weight Z-type logical

    minimise  sum(x)
    s.t.      H_X x = 0 (mod 2)                      x commutes with X checks
              L_X[j] . x = 1 (mod 2) for at least one j   x is not a stabiliser
              x binary

Mod-2 equations are linearised with integer slack variables. If the solver
finishes, the distance is exact. If it hits the time limit, its dual bound is
still a proven lower bound, which together with a witness brackets d.

For bivariate bicycle codes two symmetries make this much faster:
  * d_X = d_Z (the code is symmetric under A <-> B^T), so one side suffices;
  * translations of the torus are automorphisms, so we may assume the
    logical touches cell 0 in one of the two blocks.
"""
from __future__ import annotations

import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix, hstack, identity, vstack, lil_matrix

from .bbcode import CSSCode


def min_logical_weight(h_commute: np.ndarray, dual_logicals: np.ndarray, upper: int | None = None,
                       fix_cells: list[int] | None = None, time_limit: float = 600.0, verbose: bool = False):
    """Minimum weight of x with h_commute x = 0 and dual_logicals x != 0 (mod 2)."""
    m, n = h_commute.shape
    k = dual_logicals.shape[0]
    wrow = h_commute.sum(axis=1)
    wl = dual_logicals.sum(axis=1)
    # variable order: x (n) | s (m) | t (k) | y (k)
    nv = n + m + k + k
    H = csr_matrix(h_commute.astype(float))
    L = csr_matrix(dual_logicals.astype(float))
    A1 = hstack([H, -2 * identity(m), csr_matrix((m, 2 * k))])                      # Hx - 2s = 0
    A2 = hstack([L, csr_matrix((k, m)), -2 * identity(k), -identity(k)])             # Lx - 2t - y = 0
    A3 = hstack([csr_matrix((1, n + m + k)), csr_matrix(np.ones((1, k)))])           # sum y >= 1
    rows = [A1, A2, A3]
    lo = [np.zeros(m), np.zeros(k), np.ones(1)]
    hi = [np.zeros(m), np.zeros(k), np.full(1, np.inf)]
    if fix_cells:
        A4 = lil_matrix((1, nv))
        for q in fix_cells:
            A4[0, q] = 1
        rows.append(csr_matrix(A4))
        lo.append(np.ones(1))
        hi.append(np.full(1, np.inf))
    if upper is not None:
        A5 = hstack([csr_matrix(np.ones((1, n))), csr_matrix((1, m + 2 * k))])
        rows.append(A5)
        lo.append(np.zeros(1))
        hi.append(np.full(1, float(upper)))
    A = vstack(rows).tocsr()
    c = np.concatenate([np.ones(n), np.zeros(m + 2 * k)])
    lb = np.zeros(nv)
    ub = np.concatenate([np.ones(n), np.floor(wrow / 2), np.floor(wl / 2), np.ones(k)])
    integrality = np.ones(nv)
    t0 = time.time()
    res = milp(c, constraints=LinearConstraint(A, np.concatenate(lo), np.concatenate(hi)),
               bounds=Bounds(lb, ub), integrality=integrality,
               options={"time_limit": time_limit, "disp": verbose, "mip_rel_gap": 0})
    out = {"status": res.status, "message": res.message, "seconds": round(time.time() - t0, 1)}
    if res.x is not None:
        out["best_weight"] = int(round(res.x[:n].sum()))
        out["witness"] = np.nonzero(np.round(res.x[:n]))[0].tolist()
    lb_val = getattr(res, "mip_dual_bound", None)
    out["lower_bound"] = int(np.ceil(lb_val - 1e-6)) if lb_val is not None and np.isfinite(lb_val) else None
    out["exact"] = res.status == 0
    if out["exact"]:
        out["lower_bound"] = out["best_weight"]
    return out


def check_witness(code: CSSCode, support: list[int]) -> bool:
    """True if the support is a nontrivial Z-type logical: commutes with H_X, anticommutes with some X logical."""
    x = np.zeros(code.n, dtype=np.int64)
    x[support] = 1
    lx, _ = code.logicals
    return not ((code.hx.astype(np.int64) @ x) & 1).any() and bool(((lx.astype(np.int64) @ x) & 1).any())


def certify_bb(code: CSSCode, l: int, m: int, upper: int | None = None, time_limit: float = 600.0,
               verbose: bool = False) -> dict:
    lx, lz = code.logicals
    lm = l * m
    # Z-type logicals: commute with X checks, detected by X logicals
    return min_logical_weight(code.hx, lx, upper=upper, fix_cells=[0, lm], time_limit=time_limit, verbose=verbose)


def main(argv=None):
    """python -m qec_search.exact_distance --genome '<json>' --claimed 20 --hours 6 --out result.json"""
    import argparse
    import json
    from .bbcode import BBGenome, build_bb

    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--claimed", type=int, required=True, help="weight of a known logical (upper bound)")
    ap.add_argument("--hours", type=float, default=6.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    g = BBGenome.from_json(args.genome)
    code = build_bb(g)
    print(f"[[{code.n},{code.k},<={args.claimed}]]  {g.pretty()}")
    print("solver log follows; 'dual bound' is the proven minimum so far, it rises over time", flush=True)
    # no upper cap: scipy only reports the proven bound at the time limit if the
    # solver holds some solution, and without the cap it finds one quickly
    r = certify_bb(code, g.l, g.m, upper=None, time_limit=args.hours * 3600, verbose=True)
    # keep the witness so anyone can check the upper side without re-solving
    if "witness" in r:
        r["witness_valid"] = check_witness(code, r["witness"])
    r.update({"genome": g.to_json(), "claimed": args.claimed})
    if r["exact"]:
        verdict = f"PROVEN: minimum distance is exactly {r['best_weight']}"
    elif r.get("lower_bound"):
        verdict = f"time limit reached: distance proven >= {r['lower_bound']} (and <= {args.claimed})"
    else:
        verdict = "time limit reached with no proven bound"
    r["verdict"] = verdict
    print("\n" + verdict, flush=True)
    if args.out:
        with open(args.out, "w") as f:
            json.dump(r, f, indent=1)
        print(f"saved {args.out}")


if __name__ == "__main__":
    main()
