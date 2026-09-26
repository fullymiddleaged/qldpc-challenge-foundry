"""Reference BB codes from Bravyi et al., "High-threshold and low-overhead
fault-tolerant quantum memory", Nature 627, 778 (2024), Table 3.

These are the baselines the search is measured against. A candidate that
beats them on k*d^2/n at the same n (and holds up at circuit level) is
interesting. A candidate whose parameters are not in this list is NOT
necessarily new: many more BB and generalized-bicycle codes have been
tabulated since (Lin & Pryadko; Wang & Mueller; the 2026 LLM-guided searches).
Always check the literature before claiming novelty.
"""
from .bbcode import BBGenome

X = lambda i: (i, 0)
Y = lambda j: (0, j)

KNOWN = [
    # name, genome, (n, k, d)
    ("[[72,12,6]]", BBGenome(6, 6, (X(3), Y(1), Y(2)), (Y(3), X(1), X(2))), (72, 12, 6)),
    ("[[90,8,10]]", BBGenome(15, 3, (X(9), Y(1), Y(2)), ((0, 0), X(2), X(7))), (90, 8, 10)),
    ("[[108,8,10]]", BBGenome(9, 6, (X(3), Y(1), Y(2)), (Y(3), X(1), X(2))), (108, 8, 10)),
    ("[[144,12,12]] gross", BBGenome(12, 6, (X(3), Y(1), Y(2)), (Y(3), X(1), X(2))), (144, 12, 12)),
    ("[[288,12,18]]", BBGenome(12, 12, (X(3), Y(2), Y(7)), (Y(3), X(1), X(2))), (288, 12, 18)),
]

KNOWN_PARAMS = {p for _, _, p in KNOWN}
KNOWN_KEYS = {g.canonical_key(): name for name, g, _ in KNOWN}


def _known_sigs():
    from .bbcode import build_bb
    return {build_bb(g).signature(): name for name, g, _ in KNOWN}


KNOWN_SIGS = _known_sigs()


def best_known_kd2n(n: int) -> float:
    """Best k*d^2/n among the reference codes with at most n qubits."""
    vals = [k * d * d / nn for _, _, (nn, k, d) in KNOWN if nn <= n]
    return max(vals) if vals else 0.0
