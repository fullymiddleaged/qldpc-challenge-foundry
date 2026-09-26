"""Small GF(2) linear-algebra toolkit on numpy uint8 arrays.

Everything here is pure numpy so it runs anywhere. Matrices are 2-D arrays
of 0/1 (dtype uint8). Row operations are vectorised per pivot column, which is
fast enough for codes up to a few thousand qubits.
"""
from __future__ import annotations

import numpy as np


def as_gf2(a) -> np.ndarray:
    return (np.asarray(a) & 1).astype(np.uint8)


def rref(mat: np.ndarray, col_order: np.ndarray | None = None):
    """Row-reduce a copy of `mat` over GF(2).

    Columns are visited in `col_order` (default: left to right), which lets
    callers bias which columns become pivots (used by OSD and by the
    information-set distance search).

    Returns (reduced_matrix, pivot_columns) where reduced_matrix has the
    same shape as `mat` and only its first len(pivot_columns) rows are
    non-zero.
    """
    m = as_gf2(mat).copy()
    rows, cols = m.shape
    order = np.arange(cols) if col_order is None else np.asarray(col_order)
    pivots = []
    r = 0
    for c in order:
        if r >= rows:
            break
        nz = np.nonzero(m[r:, c])[0]
        if nz.size == 0:
            continue
        p = r + nz[0]
        if p != r:
            m[[r, p]] = m[[p, r]]
        hits = np.nonzero(m[:, c])[0]
        hits = hits[hits != r]
        if hits.size:
            m[hits] ^= m[r]
        pivots.append(int(c))
        r += 1
    return m, pivots


def rank(mat: np.ndarray) -> int:
    return len(rref(mat)[1])


def nullspace(mat: np.ndarray, col_order: np.ndarray | None = None) -> np.ndarray:
    """Basis of {v : mat @ v = 0 mod 2}, one vector per row.

    With a custom `col_order`, the columns visited first become pivots, so the
    free (information) columns are the ones visited last. Each basis vector
    has exactly one free column set, which is what the information-set
    distance search relies on.
    """
    mat = as_gf2(mat)
    n = mat.shape[1]
    red, pivots = rref(mat, col_order)
    pivot_set = set(pivots)
    free = [c for c in (range(n) if col_order is None else col_order) if c not in pivot_set]
    basis = np.zeros((len(free), n), dtype=np.uint8)
    piv_rows = red[: len(pivots)]
    for i, f in enumerate(free):
        basis[i, f] = 1
        # pivot variable in row j equals the value of column f in that row
        col = piv_rows[:, f]
        for j in np.nonzero(col)[0]:
            basis[i, pivots[j]] = 1
    return basis


def row_basis(mat: np.ndarray) -> np.ndarray:
    red, piv = rref(mat)
    return red[: len(piv)]


def solve(mat: np.ndarray, rhs: np.ndarray, col_order=None):
    """Find one x with mat @ x = rhs (mod 2), or None if inconsistent.

    Free variables are set to zero, so with a reliability-sorted `col_order`
    this is exactly OSD-0.
    """
    mat = as_gf2(mat)
    rhs = as_gf2(rhs).reshape(-1, 1)
    aug = np.hstack([mat, rhs])
    n = mat.shape[1]
    order = np.arange(n) if col_order is None else np.asarray(col_order)
    red, piv = rref(aug, order)  # never pivots on the rhs column
    rows_used = len(piv)
    if red[rows_used:, -1].any():
        return None
    x = np.zeros(n, dtype=np.uint8)
    for j, c in enumerate(piv):
        x[c] = red[j, -1]
    return x


def independent_extension(base: np.ndarray, candidates: np.ndarray) -> np.ndarray:
    """Pick rows of `candidates` that are independent of rowspace(base) and of
    each other. Used to find logical operators: ker(H_Z) modulo rowspace(H_X).
    """
    chosen = []
    current = row_basis(base) if base.size else np.zeros((0, candidates.shape[1]), np.uint8)
    r0 = current.shape[0]
    for v in candidates:
        trial = np.vstack([current, v[None, :]])
        if rank(trial) > current.shape[0]:
            current = row_basis(trial)
            chosen.append(v)
    assert current.shape[0] == r0 + len(chosen)
    return np.array(chosen, dtype=np.uint8).reshape(len(chosen), candidates.shape[1])
