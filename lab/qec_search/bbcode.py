"""Bivariate bicycle (BB) codes, the family behind IBM's gross code.

A BB code lives on an l x m torus. Let x = S_l (x) I_m and y = I_l (x) S_m,
where S_k is the k x k cyclic shift. Pick two polynomials A and B, each a sum
of monomials x^i y^j. Then

    H_X = [A | B]        H_Z = [B^T | A^T]

and because x and y commute, H_X H_Z^T = AB + BA = 0, so this is always a
valid CSS code with n = 2 l m physical qubits and check weight
|A| + |B|.

A genome is (l, m, A_terms, B_terms) with each term an (i, j) exponent pair.
IBM's published codes use three terms each, every term a pure power of x or y.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
import json

import numpy as np

from . import gf2


def _shift(k: int) -> np.ndarray:
    return np.roll(np.eye(k, dtype=np.uint8), 1, axis=1)


def _monomial(l: int, m: int, i: int, j: int) -> np.ndarray:
    xi = np.linalg.matrix_power(_shift(l).astype(np.int64), i % l) if l > 1 else np.eye(1, dtype=np.int64)
    yj = np.linalg.matrix_power(_shift(m).astype(np.int64), j % m) if m > 1 else np.eye(1, dtype=np.int64)
    return (np.kron(xi, yj) & 1).astype(np.uint8)


def poly_matrix(l: int, m: int, terms) -> np.ndarray:
    out = np.zeros((l * m, l * m), dtype=np.uint8)
    for i, j in terms:
        out ^= _monomial(l, m, i, j)
    return out


def term_str(t) -> str:
    i, j = t
    parts = []
    if i:
        parts.append("x" if i == 1 else f"x^{i}")
    if j:
        parts.append("y" if j == 1 else f"y^{j}")
    return "".join(parts) or "1"


@dataclass(frozen=True)
class BBGenome:
    l: int
    m: int
    A: tuple  # tuple of (i, j)
    B: tuple

    def normalized(self) -> "BBGenome":
        norm = lambda ts: tuple(sorted({(i % self.l, j % self.m) for i, j in ts}))
        return BBGenome(self.l, self.m, norm(self.A), norm(self.B))

    def is_valid(self) -> bool:
        g = self.normalized()
        # Terms must be distinct or they cancel mod 2 and lower the check weight.
        return len(g.A) == len(self.A) and len(g.B) == len(self.B) and self.l >= 2 and self.m >= 2

    def key(self) -> str:
        g = self.normalized()
        return f"{g.l}x{g.m}|" + ",".join(f"{i}.{j}" for i, j in g.A) + "|" + ",".join(f"{i}.{j}" for i, j in g.B)

    def canonical_key(self) -> str:
        """Key shared by codes that are trivially equivalent.

        Covers: swapping A and B; multiplying A or B by a monomial (a column
        permutation); the automorphisms x -> x^u, y -> y^v for units u, v
        (this includes inversion); and exchanging x and y (transposing the
        torus). Codes with different canonical keys can still be equivalent
        in less obvious ways, so this is a de-duplication aid, not a proof.
        """
        from math import gcd

        l, m = self.l, self.m
        def shifts(poly, L, M):
            outs = []
            for (a, b) in poly:  # make each term in turn the constant term 1
                outs.append(tuple(sorted(((i - a) % L, (j - b) % M) for i, j in poly)))
            return outs

        best = None
        variants = [(l, m, self.A, self.B)]
        sw = lambda P: tuple((j, i) for i, j in P)
        variants.append((m, l, sw(self.A), sw(self.B)))  # x <-> y exchange
        for (L, M, A, B) in variants:
            ul = [u for u in range(1, L) if gcd(u, L) == 1] or [1]
            um = [v for v in range(1, M) if gcd(v, M) == 1] or [1]
            for u in ul:
                for v in um:
                    A2 = [((i * u) % L, (j * v) % M) for i, j in A]
                    B2 = [((i * u) % L, (j * v) % M) for i, j in B]
                    for P, Q in ((A2, B2), (B2, A2)):
                        cand = (L, M, min(shifts(P, L, M)), min(shifts(Q, L, M)))
                        if best is None or cand < best:
                            best = cand
        L, M, P, Q = best
        return f"{L}x{M}|" + ",".join(f"{i}.{j}" for i, j in P) + "|" + ",".join(f"{i}.{j}" for i, j in Q)

    def pretty(self) -> str:
        g = self.normalized()
        a = " + ".join(term_str(t) for t in g.A)
        b = " + ".join(term_str(t) for t in g.B)
        return f"l={g.l}, m={g.m}, A = {a}, B = {b}"

    def to_json(self) -> dict:
        g = self.normalized()
        return {"l": g.l, "m": g.m, "A": [list(t) for t in g.A], "B": [list(t) for t in g.B]}

    @staticmethod
    def from_json(d) -> "BBGenome":
        if isinstance(d, str):
            d = json.loads(d)
        return BBGenome(int(d["l"]), int(d["m"]), tuple(tuple(t) for t in d["A"]), tuple(tuple(t) for t in d["B"]))


@dataclass
class CSSCode:
    hx: np.ndarray
    hz: np.ndarray
    name: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return self.hx.shape[1]

    @cached_property
    def rank_x(self) -> int:
        return gf2.rank(self.hx)

    @cached_property
    def rank_z(self) -> int:
        return gf2.rank(self.hz)

    @property
    def k(self) -> int:
        return self.n - self.rank_x - self.rank_z

    def check_css(self) -> bool:
        return not ((self.hx.astype(np.int64) @ self.hz.T.astype(np.int64)) & 1).any()

    @cached_property
    def logicals(self):
        """(LX, LZ): X-type logicals in ker(H_Z) / rowspace(H_X), Z-type in
        ker(H_X) / rowspace(H_Z). Each has k rows."""
        lx = gf2.independent_extension(self.hx, gf2.nullspace(self.hz))
        lz = gf2.independent_extension(self.hz, gf2.nullspace(self.hx))
        return lx, lz

    def components(self) -> int:
        """Number of connected pieces of the Tanner graph. More than 1 means
        the code is several independent smaller codes side by side (e.g. two
        copies of the gross code), which is never interesting."""
        import scipy.sparse as sp
        from scipy.sparse.csgraph import connected_components
        h = sp.csr_matrix(np.vstack([self.hx, self.hz]))
        adj = sp.bmat([[None, h], [h.T, None]])
        return int(connected_components(adj, directed=False)[0])

    def signature(self) -> str:
        """Fingerprint that is identical for codes that are the same up to
        relabelling qubits and checks (and swapping X/Z). Built from the
        eigenvalues of H^T H. Different fingerprints prove two codes are
        different; equal fingerprints almost always mean the same code."""
        import hashlib
        parts = []
        for h in (self.hx, self.hz):
            ev = np.round(np.linalg.eigvalsh(h.T.astype(float) @ h.astype(float)), 6) + 0.0
            parts.append(hashlib.md5(ev.tobytes()).hexdigest()[:12])
        return f"{self.n}.{self.k}." + "-".join(sorted(parts))

    def distance_upper_bound(self, trials: int = 300, seed: int = 0, which: str = "both",
                             stop_at: int | None = None) -> int:
        """Randomised information-set search for low-weight logicals.

        Returns the lowest-weight non-trivial logical found, which is an upper
        bound on the true distance (it equals it with high probability once
        trials are large enough, but it is not a proof).
        """
        rng = np.random.default_rng(seed)
        lx, lz = self.logicals
        best = self.n
        targets = []
        if which in ("both", "z"):
            targets.append((self.hx, lx))  # Z-type logicals: in ker(H_X), detected via LX
        if which in ("both", "x"):
            targets.append((self.hz, lz))
        for _ in range(trials):
            for h, dual in targets:
                perm = rng.permutation(self.n)
                basis = gf2.nullspace(h, perm)
                w = basis.sum(axis=1)
                nontriv = ((basis.astype(np.int64) @ dual.T.astype(np.int64)) & 1).any(axis=1)
                if nontriv.any():
                    best = min(best, int(w[nontriv].min()))
                # pairwise sums of the lightest basis vectors catch a few more
                light = np.argsort(w)[:12]
                for a in range(len(light)):
                    for b in range(a + 1, len(light)):
                        v = basis[light[a]] ^ basis[light[b]]
                        if ((dual.astype(np.int64) @ v) & 1).any():
                            best = min(best, int(v.sum()))
            if stop_at is not None and best <= stop_at:
                break
        return best


def build_bb(g: BBGenome) -> CSSCode:
    g = g.normalized()
    a = poly_matrix(g.l, g.m, g.A)
    b = poly_matrix(g.l, g.m, g.B)
    hx = np.hstack([a, b])
    hz = np.hstack([b.T, a.T])
    return CSSCode(hx=hx, hz=hz, name=g.pretty(), meta={"genome": g.to_json()})
