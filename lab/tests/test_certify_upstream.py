from __future__ import annotations

import importlib.util
import itertools
from pathlib import Path

import numpy as np

from qec_search.bbcode import BBGenome, build_bb

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("certify_upstream", ROOT / "scripts/certify_upstream.py")
cu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cu)

GROSS = BBGenome(l=12, m=6, A=((3, 0), (0, 1), (0, 2)), B=((0, 3), (1, 0), (2, 0)))


def _matches(cube, v):
    return all(v[abs(x) - 1] == (x > 0) for x in cube)


def test_torus_shape():
    assert cu.torus_shape(Path("bb_216_8_16_18x6.npz")) == (18, 6)


def test_root_cubes_partition_the_vectors_touching_cell_zero():
    l, m = 2, 2                                          # n = 8, small enough to enumerate
    roots = cu.root_cubes(l, m, bits=2)
    assert len(roots) == 8
    for v in itertools.product((0, 1), repeat=2 * l * m):
        hits = sum(_matches(c, v) for c in roots)
        assert hits == (1 if v[0] or v[l * m] else 0)


def test_children_split_on_the_lowest_free_qubit():
    assert cu.children([1, -3], 8) == [[1, -3, 2], [1, -3, -2]]


def test_leaves_and_summary_follow_splits():
    roots = cu.root_cubes(2, 2, bits=0)                  # [[1], [-1, 5]]
    st = {"sides": ["X"], "cubes": {
        cu.cube_key("X", [1]): {"status": "SPLIT", "secs": 10},
        cu.cube_key("X", [1, 2]): {"status": "UNSAT", "secs": 1},
        cu.cube_key("X", [-1, 5]): {"status": "UNSAT", "secs": 1},
    }}
    assert [c for _, c in cu.pending(st, roots, 8)] == [[1, -2]]
    s = cu.summarize(st, roots, 8)
    assert s["verdict"] == "in progress" and s["proven_share"] == 0.75 and s["splits"] == 1
    st["cubes"][cu.cube_key("X", [1, -2])] = {"status": "UNSAT", "secs": 1}
    assert cu.summarize(st, roots, 8)["verdict"] == "PROVEN" and cu.pending(st, roots, 8) == []
    st["cubes"][cu.cube_key("X", [1, -2])] = {"status": "SAT", "secs": 1, "witness": [0]}
    assert cu.summarize(st, roots, 8)["verdict"] == "REFUTED" and cu.pending(st, roots, 8) == []


def test_gross_code_has_both_symmetries_and_a_broken_code_does_not():
    code = build_bb(GROSS)
    hx, hz = code.hx.astype(np.int8), code.hz.astype(np.int8)
    assert cu.has_translation_symmetry(hx, hz, 12, 6)
    assert cu.has_xz_duality(hx, hz, 12, 6)
    shuffled = hx[:, np.random.default_rng(0).permutation(hx.shape[1])]
    assert not cu.has_translation_symmetry(shuffled, hz, 12, 6)
