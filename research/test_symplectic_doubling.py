"""Round-trip tests for symplectic doubling and its inverse (halving).

Companion test for fieldnotes/2026-09-28-symplectic-doubling-and-halving.md.

The fieldnote claims:

  * doubling any isotropic stabilizer code S = (A | B) as H'_X = (A | B),
    H'_Z = (B | A) is CSS, with n -> 2n and k -> 2k (so kd^2/n is non-decreasing:
    d' >= d, with equality iff a minimum-weight logical is Y-free);
  * the inverse (halving) is conditional: it needs the CPM pair-partition
    structure, whose off-diagonal difference arrays must take every value an
    even number of times;
  * the fold is genuinely non-CSS yet still doubles to a CSS code.

These are cheap algebraic claims, so the test checks them directly on a small
stabilizer code (the [[5,1,3]] perfect code, for the unconditional direction)
and on the published (4,10,29) instance whose fold and parent shipped as
[[145,32,6]] / [[290,64,6]]. The builder under test is research/build_halved_pp.py.

No distance search runs here: distances in this construction are witness-backed
upper bounds, not certified here.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.join(HERE, "kit"), HERE, os.path.join(os.path.dirname(HERE), "verify")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import build_halved_pp as bh  # noqa: E402
from css import compute_k, rank, verify_css  # noqa: E402
from heuristic_distance import doubled_matrices  # noqa: E402
from qldpc_verify import is_css_up_to_local_hadamard  # noqa: E402

# Published instance: (J, L, P) = (4, 10, 29), sigma(l) = (l + 5) mod 10,
# eta = -1. Its parent [[290,64,6]] and fold [[145,32,6]] are the worked pair.
P_PUB, J_PUB, L_PUB, ETA_PUB = 29, 4, 10, -1
E_PUB = np.array(
    [
        [2, 20, 21, 27, 6, 21, 13, 12, 6, 14],
        [11, 22, 23, 0, 15, 23, 12, 8, 8, 19],
        [5, 26, 8, 1, 3, 8, 26, 5, 22, 1],
        [12, 23, 24, 24, 16, 21, 3, 12, 25, 17],
    ]
)


def _published():
    """Build the published CSS parent, its D array, and the shift involution."""
    sig = bh.make_sig("shift", L_PUB)
    D = np.array([[ETA_PUB * int(E_PUB[j, sig(ell)]) % P_PUB for ell in range(L_PUB)] for j in range(J_PUB)])
    HX, HZ = bh.circulant_build(E_PUB, D, P_PUB)
    return HX, HZ, D, sig


def _perfect_5_1_3():
    """[[5,1,3]] as S = (A | B): the four cyclic XZZXI stabilizers."""
    A = np.array([[1, 0, 0, 1, 0], [0, 1, 0, 0, 1], [1, 0, 1, 0, 0], [0, 1, 0, 1, 0]], dtype=np.int8)
    B = np.array([[0, 1, 1, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 1, 1], [1, 0, 0, 0, 1]], dtype=np.int8)
    return A, B


def test_doubling_is_always_css_and_doubles_k():
    """The forward map has no precondition: any isotropic S doubles to CSS."""
    A, B = _perfect_5_1_3()
    assert not ((A @ B.T + B @ A.T) % 2).any(), "[[5,1,3]] must be isotropic"
    n = A.shape[1]
    k = n - rank(np.concatenate([A, B], axis=1))
    assert k == 1

    HX2, HZ2 = doubled_matrices(A, B)
    assert verify_css(HX2, HZ2)  # CSS commutation holds
    assert HX2.shape == (A.shape[0], 2 * n)  # n -> 2n
    assert compute_k(HX2, HZ2) == 2 * k  # k -> 2k


def test_halving_conditions_hold_for_the_published_instance():
    """Check the inverse's conditions on the published instance.

    D = eta * E∘sigma holds, and every off-diagonal difference array takes each
    value an even number of times.
    """
    _, _, D, sig = _published()
    for j in range(J_PUB):
        for ell in range(L_PUB):
            assert D[j, ell] == (ETA_PUB * int(E_PUB[j, sig(ell)])) % P_PUB
    for i in range(J_PUB):
        for j in range(J_PUB):
            if i == j:
                continue
            counts = {}
            for ell in range(L_PUB):
                v = (int(E_PUB[i, ell]) - ETA_PUB * int(E_PUB[j, sig(ell)])) % P_PUB
                counts[v] = counts.get(v, 0) + 1
            assert all(c % 2 == 0 for c in counts.values()), (i, j, counts)


def test_published_pair_round_trips():
    """Parent -> fold (halve) -> double reproduces a CSS code with the parent k."""
    HX, HZ, _, sig = _published()
    assert verify_css(HX, HZ)
    k_parent = compute_k(HX, HZ)
    assert (HX.shape[1], k_parent) == (L_PUB * P_PUB, 64)

    S, n_fold, k_fold = bh.halve(HX, HZ, J_PUB, L_PUB, P_PUB, eta=ETA_PUB, sig=sig)
    assert n_fold == L_PUB * P_PUB // 2 == 145
    assert k_fold == k_parent // 2 == 32

    A, B = S[:, :n_fold], S[:, n_fold:]
    assert not ((A @ B.T + B @ A.T) % 2).any(), "the fold must be isotropic"

    HX2, HZ2 = doubled_matrices(A, B)
    assert verify_css(HX2, HZ2)
    assert HX2.shape[1] == 2 * n_fold == 290
    assert compute_k(HX2, HZ2) == k_parent == 64


def test_fold_is_not_css_up_to_local_hadamard():
    """The fold is irreducibly non-CSS, yet it still doubles to a CSS code."""
    HX, HZ, _, sig = _published()
    S, n_fold, _ = bh.halve(HX, HZ, J_PUB, L_PUB, P_PUB, eta=ETA_PUB, sig=sig)
    A, B = S[:, :n_fold], S[:, n_fold:]
    assert is_css_up_to_local_hadamard(A, B) is None
    assert verify_css(*doubled_matrices(A, B))


def _combined_components(HX, HZ, n):
    """Count connected components of the combined X/Z Tanner graph."""
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for M in (HX, HZ):
        for row in M:
            qubits = [int(q) for q in np.nonzero(row)[0]]
            for q in qubits[1:]:
                ra, rb = find(qubits[0]), find(q)
                if ra != rb:
                    parent[ra] = rb
    return len({find(i) for i in range(n)})


def test_doubling_a_css_input_is_a_disconnected_direct_sum():
    """CSS in, CSS out is degenerate: the double is two disconnected copies."""
    # Steane [[7,1,3]] in canonical pure form: A = [H; 0], B = [0; H].
    H = np.array([[1, 1, 1, 0, 1, 0, 0], [1, 0, 1, 1, 0, 1, 0], [0, 1, 1, 1, 0, 0, 1]], dtype=np.int8)
    A = np.vstack([H, np.zeros_like(H)])
    B = np.vstack([np.zeros_like(H), H])
    HX2, HZ2 = doubled_matrices(A, B)
    assert verify_css(HX2, HZ2)
    assert compute_k(HX2, HZ2) == 2
    assert _combined_components(HX2, HZ2, 14) == 2  # direct sum, inadmissible

    # Contrast: a genuinely non-CSS input doubles to a connected CSS code.
    A5, B5 = _perfect_5_1_3()
    HX5, HZ5 = doubled_matrices(A5, B5)
    assert verify_css(HX5, HZ5)
    assert _combined_components(HX5, HZ5, 10) == 1
