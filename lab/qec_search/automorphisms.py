# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Automorphisms of a CSS code's Tanner graph, found by individualisation-refinement.

An automorphism here is a qubit permutation p that maps the set of X-check supports onto itself and the set of
Z-check supports onto itself. A duality is a qubit permutation that maps the X-check supports onto the Z-check
supports and vice versa. Either kind preserves the relevant rowspaces, so it maps minimum-weight logicals to
minimum-weight logicals (a duality swaps the X and Z sides, which proves d_X = d_Z). These are the Tanner-graph
automorphisms, a subgroup of the code's full automorphism group. That is enough for symmetry breaking, since any
subgroup is sound to branch on.

The search is the usual one behind nauty-style tools, kept small. Both graphs go into one disjoint union and are
colour-refined (1-WL) together, so a colour means the same thing on both sides. A vertex of the smallest
non-singleton cell on the left is paired with each same-coloured vertex on the right, and the search recurses,
backtracking on failure. The search is complete: `find_isomorphism` returns a mapping whenever one exists with the
requested pairs. Every returned permutation is also checked directly against the check supports, so a refinement bug
could only make the search miss symmetries, never report a false one.

Permutations are integer arrays p with qubit q -> p[q].
"""
from __future__ import annotations

import numpy as np

_QUBIT, _X, _Z = 0, 1, 2


def _rows(h: np.ndarray) -> list[tuple[int, ...]]:
    return [tuple(int(j) for j in np.nonzero(r)[0]) for r in np.asarray(h) if r.any()]


class TannerGraph:
    """Qubits 0..n-1, then X checks, then Z checks, as adjacency lists with an initial vertex colour."""

    def __init__(self, hx: np.ndarray, hz: np.ndarray, swap: bool = False):
        self.n = int(np.asarray(hx).shape[1])
        self.x_rows, self.z_rows = _rows(hx), _rows(hz)
        adj: list[list[int]] = [[] for _ in range(self.n)]
        colour = [_QUBIT] * self.n
        for kind, rows in ((_X, self.x_rows), (_Z, self.z_rows)):
            for row in rows:
                v = len(adj)
                adj.append(list(row))
                for q in row:
                    adj[q].append(v)
                colour.append({_X: _Z, _Z: _X}[kind] if swap else kind)
        self.adj, self.colour = adj, colour


def _refine(adj: list[list[int]], colour: list[int]) -> list[int]:
    """Colour refinement to the coarsest equitable partition; returns canonical colour ids."""
    cur = colour
    ncls = len(set(cur))
    while True:
        sig = [(cur[v], tuple(sorted(cur[u] for u in adj[v]))) for v in range(len(adj))]
        ids = {s: i for i, s in enumerate(sorted(set(sig)))}
        nxt = [ids[s] for s in sig]
        if len(ids) == ncls:
            return nxt
        cur, ncls = nxt, len(ids)


def _union(a: TannerGraph, b: TannerGraph) -> tuple[list[list[int]], list[int], int]:
    off = len(a.adj)
    adj = [list(nb) for nb in a.adj] + [[u + off for u in nb] for nb in b.adj]
    return adj, list(a.colour) + list(b.colour), off


def find_isomorphism(a: TannerGraph, b: TannerGraph, pairs: list[tuple[int, int]] = (),
                     budget: int = 2000) -> np.ndarray | None:
    """A qubit permutation mapping graph a onto b and each qubit pairs[i][0] to pairs[i][1], or None.

    `budget` caps the refinements tried; running out returns None. That can only lose symmetries, and branching on
    a smaller group is still sound."""
    if a.n != b.n or len(a.adj) != len(b.adj):
        return None
    adj, base, off = _union(a, b)
    total = len(adj)
    left_budget = [budget]

    def search(colour: list[int]) -> list[int] | None:
        if left_budget[0] <= 0:
            return None
        left_budget[0] -= 1
        col = _refine(adj, colour)
        left: dict[int, list[int]] = {}
        right: dict[int, list[int]] = {}
        for v in range(total):
            (left if v < off else right).setdefault(col[v], []).append(v)
        if left.keys() != right.keys() or any(len(left[c]) != len(right[c]) for c in left):
            return None
        cells = [c for c in left if len(left[c]) > 1]
        if not cells:
            return [right[col[v]][0] - off for v in range(off)]
        c = min(cells, key=lambda c: len(left[c]))
        u = left[c][0]
        fresh = max(col) + 1
        for w in right[c]:
            trial = list(col)
            trial[u] = trial[w] = fresh
            found = search(trial)
            if found is not None:
                return found
        return None

    colour = list(base)
    fresh = max(colour) + 1
    for i, (p, q) in enumerate(pairs):
        colour[p] = colour[q + off] = fresh + i
    full = search(colour)
    if full is None:
        return None
    perm = np.array(full[: a.n], dtype=np.int64)
    return perm if maps_checks(a, b, perm) else None


def maps_checks(a: TannerGraph, b: TannerGraph, perm: np.ndarray) -> bool:
    """Direct check: perm sends a's X supports onto b's X supports and a's Z supports onto b's Z supports.

    For a graph built with swap=True, b's "X" colour marks the original Z checks, which is what makes a duality."""
    b_by_colour = {_X: [], _Z: []}
    rows_b = b.x_rows + b.z_rows
    for i, row in enumerate(rows_b):
        b_by_colour[b.colour[b.n + i]].append(tuple(sorted(row)))
    rows_a = a.x_rows + a.z_rows
    a_by_colour = {_X: [], _Z: []}
    for i, row in enumerate(rows_a):
        a_by_colour[a.colour[a.n + i]].append(tuple(sorted(int(perm[q]) for q in row)))
    return all(sorted(a_by_colour[k]) == sorted(b_by_colour[k]) for k in (_X, _Z))


def _orbits_from_generators(n: int, gens: list[np.ndarray]) -> list[list[int]]:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in gens:
        for x in range(n):
            rx, ry = find(x), find(int(g[x]))
            if rx != ry:
                parent[max(rx, ry)] = min(rx, ry)
    out: dict[int, list[int]] = {}
    for x in range(n):
        out.setdefault(find(x), []).append(x)
    return sorted(out.values(), key=lambda o: o[0])


def qubit_orbits(hx: np.ndarray, hz: np.ndarray, fixed: tuple[int, ...] = (),
                 among: set[int] | None = None) -> tuple[list[list[int]], list[np.ndarray]]:
    """Orbits of the qubits (or of `among`) under the automorphisms that fix every qubit in `fixed` pointwise, with
    the generators found.

    With `among`, a generator is kept only if it maps `among` onto itself, so the group generated preserves that set
    and its complement. Orbits are exact for the group generated. It equals the whole pointwise stabiliser unless the
    search budget ran out, and a subgroup is still sound to branch on."""
    g = TannerGraph(hx, hz)
    n = g.n
    base = list(g.colour)
    fresh = max(base) + 1
    for i, q in enumerate(fixed):
        base[q] = fresh + i
    cell = _refine(g.adj, base)
    pool = sorted(among) if among is not None else [q for q in range(n) if q not in fixed]
    keep = set(pool)
    pool_arr = np.array(pool, dtype=np.int64)
    gens: list[np.ndarray] = []
    for q in pool:
        orbits = _orbits_from_generators(n, gens)
        of = {x: min(o) for o in orbits for x in o}
        reps = sorted({of[r] for r in pool if r < q and cell[r] == cell[q]})
        if of[q] in reps:
            continue
        for r in reps:
            perm = find_isomorphism(g, g, [(f, f) for f in fixed] + [(r, q)])
            if perm is not None and (among is None or set(perm[pool_arr].tolist()) == keep):
                gens.append(perm)
                break
    orbits = [[x for x in o if x in keep] for o in _orbits_from_generators(n, gens)]
    return [o for o in orbits if o], gens


def find_duality(hx: np.ndarray, hz: np.ndarray) -> np.ndarray | None:
    """A qubit permutation swapping the X and Z check sets, or None."""
    return find_isomorphism(TannerGraph(hx, hz), TannerGraph(hx, hz, swap=True))
