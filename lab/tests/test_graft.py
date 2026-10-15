# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
from __future__ import annotations

import numpy as np
import pytest

from qec_search import graft, tile
from qec_search.bbcode import CSSCode

A637 = [("h", 0, 0), ("h", 0, 1), ("h", 1, 1), ("h", 2, 2), ("v", 0, 2), ("v", 1, 0), ("v", 2, 0), ("v", 2, 1)]


@pytest.fixture(scope="module")
def small_tile():
    hx, hz, qubits = tile.build(A637, 3, 5, 5)
    return hx, hz, tile.layout(qubits)


def test_low_weight_elements_are_row_sums_lightest_first(small_tile):
    hx, _, _ = small_tile
    elems = graft.low_weight_elements(hx, max_terms=3, max_weight=8)
    sizes = [len(s) for s, _ in elems]
    assert sizes == sorted(sizes) and sizes[0] >= 1
    for support, subset in elems[:200]:
        assert 1 <= len(subset) <= 3
        assert np.array_equal(np.nonzero(hx[list(subset)].sum(0) % 2)[0], np.array(support))
    assert len({s for s, _ in elems}) == len(elems)


def test_apply_deletes_one_qubit_and_keeps_the_code_valid(small_tile):
    hx, hz, _ = small_tile
    k = graft.k_of(hx, hz)
    support, subset = graft.low_weight_elements(hx, max_terms=2)[0]
    a = support[0]
    nx, nz, keep, _ = graft.apply(hx, hz, support, subset, a)
    assert nx.shape == (hx.shape[0] - 1, hx.shape[1] - 1) and nz.shape == (hz.shape[0], hz.shape[1] - 1)
    assert a not in keep
    assert not ((nx.astype(int) @ nz.T.astype(int)) % 2).any()
    assert graft.k_of(nx, nz) == k
    assert (nz.sum(1) <= hz.sum(1)).all()                      # Z rows only lose an index


def test_apply_refuses_a_support_outside_the_row_space(small_tile):
    hx, hz, _ = small_tile
    bad = next(q for q in range(hx.shape[1]) if hz[:, q].any())
    with pytest.raises(ValueError):
        graft.apply(hx, hz, (bad,), (0,), bad)


def test_reduce_strips_a_weight_one_stabilizer_without_a_distance_check(small_tile):
    hx, hz, coords = small_tile
    n = hx.shape[1]
    px = np.zeros((hx.shape[0] + 1, n + 1), dtype=np.uint8)
    px[:-1, :n] = hx
    px[-1, n] = 1                                              # extra qubit fixed by X_n alone
    pz = np.hstack([hz, np.zeros((hz.shape[0], 1), dtype=np.uint8)])
    pc = np.vstack([coords, coords.max(0) + 1])
    calls = []
    rx, rz, _, idx = graft.reduce(px, pz, pc, w_cap=8, r_cap=4.0, dist_ok=lambda *_: calls.append(1) or False,
                                  log=lambda *_: None)
    assert list(idx) == list(range(n)) and calls                # only the free graft went through
    assert np.array_equal(rx, hx) and np.array_equal(rz, hz)


def test_reduce_keeps_k_weight_radius_and_distance(small_tile):
    hx, hz, coords = small_tile
    code = CSSCode(hx=hx, hz=hz, name="t", meta={})
    d = code.distance_upper_bound(trials=200)

    def dist_ok(nx, nz):
        return CSSCode(hx=nx, hz=nz, name="g", meta={}).distance_upper_bound(trials=200, stop_at=d - 1) >= d

    rx, rz, rc, idx = graft.reduce(hx, hz, coords, w_cap=8, r_cap=4.0, dist_ok=dist_ok, log=lambda *_: None)
    assert graft.k_of(rx, rz) == code.k and len(idx) == rx.shape[1] < hx.shape[1]
    assert max(rx.sum(1).max(), rz.sum(1).max()) <= 8
    assert max(graft.radius(r, rc) for r in np.vstack([rx, rz])) <= 4.0 + 1e-9
    assert CSSCode(hx=rx, hz=rz, name="g", meta={}).distance_upper_bound(trials=500) >= d


def test_last_good_finds_the_last_passing_state_by_bisection():
    for n_states in range(1, 12):
        for first_bad in range(1, n_states + 1):
            seen = []
            ok = lambda i: seen.append(i) or i < first_bad          # noqa: E731
            assert graft.last_good(n_states, ok) == min(first_bad, n_states) - 1
            assert len(seen) <= max(1, n_states.bit_length())


def test_reduce_reports_every_accepted_state(small_tile):
    hx, hz, coords = small_tile
    states = []
    rx, _, _, idx = graft.reduce(hx, hz, coords, w_cap=8, r_cap=4.0, dist_ok=lambda *_: True, log=lambda *_: None,
                                 on_step=lambda *st: states.append(st))
    assert len(states) == hx.shape[1] - rx.shape[1] and np.array_equal(states[-1][3], idx)
