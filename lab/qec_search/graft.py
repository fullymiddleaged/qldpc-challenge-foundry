# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Row-space graft: delete a qubit from a CSS code at unchanged k (@dorakingx's move, upstream notes/280-8-15.md).

Take S in the row space of H_X with support T and a qubit a in T. Replace one generator of the subset summing to S by
S (the span is unchanged), add S into every other X row meeting a, and fan CNOTs out from a to T \\ {a}. S becomes
the weight-1 stabilizer X_a, column a of H_Z clears (each Z row meets T evenly), and qubit a is deleted: n -> n - 1,
k unchanged. Z rows only lose an index; an X row R meeting a becomes R + S minus a, which can grow when |S| >= 3, so
the touched rows are checked against the weight and radius caps. |S| = 1 keeps the distance; otherwise the caller's
distance check decides. The Z-side move is the same on (H_Z, H_X).

Pure numpy, so it runs in WSL beside the gf2_fast RIS engine (scripts/ris_ladder.py).
"""
from __future__ import annotations

from collections.abc import Callable, Iterator

import numpy as np

from . import gf2


def k_of(hx: np.ndarray, hz: np.ndarray) -> int:
    return hx.shape[1] - gf2.rank(hx) - gf2.rank(hz)


def low_weight_elements(h: np.ndarray, max_terms: int = 3, max_weight: int | None = None):
    """Distinct nonzero sums of <= max_terms rows of h whose overlap graph is connected, lightest first, as
    (support tuple, row subset)."""
    rows = [frozenset(np.nonzero(r)[0].tolist()) for r in h]
    on: dict[int, set[int]] = {}
    for i, r in enumerate(rows):
        for q in r:
            on.setdefault(q, set()).add(i)
    nbr = [set().union(*(on[q] for q in r)) - {i} for i, r in enumerate(rows)]
    seen: dict[frozenset, tuple[int, ...]] = {}

    def add(sub: tuple[int, ...]) -> None:
        s: frozenset = frozenset()
        for i in sub:
            s = s ^ rows[i]
        if s and (max_weight is None or len(s) <= max_weight) and s not in seen:
            seen[s] = sub

    frontier = {(i,) for i in range(len(rows)) if rows[i]}
    for sub in frontier:
        add(sub)
    for _ in range(max_terms - 1):
        grown = set()
        for sub in frontier:
            for j in set().union(*(nbr[i] for i in sub)) - set(sub):
                grown.add(tuple(sorted((*sub, j))))
        for sub in grown:
            add(sub)
        frontier = grown
    return sorted(((tuple(sorted(s)), sub) for s, sub in seen.items()), key=lambda e: (len(e[0]), e[0]))


def apply(ha: np.ndarray, hb: np.ndarray, support, subset, a: int):
    """The graft with S in the row space of ha. Returns (ha', hb', kept columns, rows of ha' that changed)."""
    s = np.zeros(ha.shape[1], dtype=np.uint8)
    s[list(support)] = 1
    if (hb[:, list(support)].sum(1) % 2).any():
        raise ValueError("S is not in the row space of ha: it fails to commute with hb")
    ha = ha.copy()
    g = subset[0]
    ha[g] = s
    meet = [int(r) for r in np.nonzero(ha[:, a])[0] if r != g]
    ha[meet] ^= s
    keep = np.array([c for c in range(ha.shape[1]) if c != a])
    rows = np.array([r for r in range(ha.shape[0]) if r != g])
    touched = [int(np.searchsorted(rows, r)) for r in meet]
    return ha[np.ix_(rows, keep)], hb[:, keep], keep, touched


def radius(row: np.ndarray, coords: np.ndarray) -> float:
    """Largest pairwise distance between the row's qubits (the board's interaction radius for one check)."""
    p = coords[np.nonzero(row)[0]]
    if len(p) < 2:
        return 0.0
    return float(np.sqrt(((p[:, None] - p[None]) ** 2).sum(-1)).max())


def candidates(hx: np.ndarray, hz: np.ndarray, coords: np.ndarray, w_cap: int, r_cap: float,
               max_terms: int = 3) -> Iterator[tuple[str, np.ndarray, np.ndarray, np.ndarray, int]]:
    """Every graft (side, hx', hz', kept columns, |S|) that keeps k, the weight cap and the radius cap, lightest S
    first, X side before Z side at each weight."""
    k = k_of(hx, hz)
    elems = {"X": low_weight_elements(hx, max_terms, w_cap), "Z": low_weight_elements(hz, max_terms, w_cap)}
    for size in sorted({len(s) for e in elems.values() for s, _ in e}):
        for side in ("X", "Z"):
            ha, hb = (hx, hz) if side == "X" else (hz, hx)
            for support, subset in (e for e in elems[side] if len(e[0]) == size):
                for a in support:
                    na, nb, keep, touched = apply(ha, hb, support, subset, a)
                    c = coords[keep]
                    if any(na[r].sum() > w_cap or radius(na[r], c) > r_cap + 1e-9 for r in touched):
                        continue
                    nx, nz = (na, nb) if side == "X" else (nb, na)
                    if k_of(nx, nz) != k:
                        continue
                    yield side, nx, nz, keep, size


def reduce(hx: np.ndarray, hz: np.ndarray, coords: np.ndarray, w_cap: int, r_cap: float,
           dist_ok: Callable[[np.ndarray, np.ndarray], bool], max_terms: int = 3, log=print):
    """Apply grafts to a fixpoint, taking the first that passes at each step. Returns (hx, hz, coords, original
    indices of the surviving qubits)."""
    idx = np.arange(hx.shape[1])
    while True:
        for side, nx, nz, keep, size in candidates(hx, hz, coords, w_cap, r_cap, max_terms):
            if size == 1 or dist_ok(nx, nz):
                gone = int(idx[np.setdiff1d(np.arange(len(idx)), keep)][0])
                log(f"graft {side} |S|={size}: qubit {gone} removed, n {hx.shape[1]} -> {nx.shape[1]}")
                hx, hz, coords, idx = nx, nz, coords[keep], idx[keep]
                break
        else:
            return hx, hz, coords, idx
