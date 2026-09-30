# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Tile codes (Steffan, Choe, Breuckmann, Pereira, Eberhardt, arXiv:2504.09171): open-boundary planar CSS codes.

A tile is a set of edges of a B x B box of cells, written ("h" | "v", x, y) with 0 <= x, y <= B - 1: h(x, y) joins
vertices (x, y) and (x + 1, y), v(x, y) joins (x, y) and (x, y + 1). The X-tile determines the Z-tile (condition T2):
h(x, y) in the X-tile gives v(B-1-x, B-1-y) in the Z-tile and vice versa, so X and Z tiles overlap evenly anywhere.

On an l x m bulk the code is built as the paper describes: both tiles at every bulk position; the physical qubits are
every edge of every bulk box, 2 (l + B - 1)(m + B - 1) of them; B - 1 rows of X-only stabilizers below and above the
bulk and B - 1 columns of Z-only stabilizers left and right, truncated to the physical qubits; empty stabilizers
dropped. Checked against the board: the tile {h00,h03,h22,h30,v01,v11,v20,v33} on 14 x 14 reproduces codes/578-18-20
(n = 578, k = 18, same spectral fingerprint).

layout() places each edge on its own integer site of the 45-degree rotated lattice (edge midpoints scaled by sqrt 2),
one qubit per site at spacing >= 1.
"""
from __future__ import annotations

import numpy as np

Edge = tuple[str, int, int]


def z_tile(x_tile: list[Edge], B: int) -> list[Edge]:
    return [("v" if k == "h" else "h", B - 1 - x, B - 1 - y) for k, x, y in x_tile]


def _placed(tile: list[Edge], px: int, py: int) -> set[Edge]:
    return {(k, px + x, py + y) for k, x, y in tile}


def build(x_tile: list[Edge], B: int, l: int, m: int) -> tuple[np.ndarray, np.ndarray, list[Edge]]:
    """(hx, hz, qubits) of the open-boundary tile code on an l x m bulk; qubits[i] is the edge of column i."""
    zt = z_tile(x_tile, B)
    bulk = [(px, py) for px in range(l) for py in range(m)]
    qubits = sorted({(k, px + a, py + b) for px, py in bulk for k in "hv" for a in range(B) for b in range(B)})
    q = set(qubits)
    xs = [_placed(x_tile, *p) for p in bulk] + [
        _placed(x_tile, px, py) & q for px in range(l) for py in [*range(-(B - 1), 0), *range(m, m + B - 1)]]
    zs = [_placed(zt, *p) for p in bulk] + [
        _placed(zt, px, py) & q for py in range(m) for px in [*range(-(B - 1), 0), *range(l, l + B - 1)]]
    idx = {e: i for i, e in enumerate(qubits)}

    def mat(sets):
        sets = [s for s in sets if s]
        h = np.zeros((len(sets), len(qubits)), dtype=np.uint8)
        for r, s in enumerate(sets):
            h[r, [idx[e] for e in s]] = 1
        return h
    return mat(xs), mat(zs), qubits


def layout(qubits: list[Edge]) -> np.ndarray:
    """Integer site per edge: midpoints (x + 1/2, y) and (x, y + 1/2) rotated by 45 degrees and scaled by sqrt 2."""
    pts = []
    for k, x, y in qubits:
        mx, my = (x + 0.5, y) if k == "h" else (x, y + 0.5)
        pts.append((mx + my - 0.5, mx - my + 0.5))
    c = np.array(pts, dtype=float)
    return (c - c.min(axis=0)).round().astype(int)


def valid_tile(x_tile: list[Edge], B: int) -> bool:
    """Distinct edges inside the box (condition T1)."""
    return len(set(x_tile)) == len(x_tile) and all(k in "hv" and 0 <= x < B and 0 <= y < B for k, x, y in x_tile)
