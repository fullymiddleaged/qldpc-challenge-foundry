# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""qec_search.tile against the board: the B = 4 tile of codes/578-18-20 on a 14 x 14 bulk must reproduce that code."""
from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

from qec_search import gf2, tile
from qec_search.bbcode import CSSCode
from qec_search.layout import check_supports, validate

BOARD = pathlib.Path(__file__).resolve().parents[2] / ".worktrees" / "frontier" / "codes"
T578 = [("h", 0, 0), ("h", 0, 3), ("h", 2, 2), ("h", 3, 0), ("v", 0, 1), ("v", 1, 1), ("v", 2, 0), ("v", 3, 3)]


def test_z_tile_overlaps_every_x_tile_evenly():
    zt = tile.z_tile(T578, 4)
    for dx in range(-4, 5):
        for dy in range(-4, 5):
            xs = {(k, x + dx, y + dy) for k, x, y in T578}
            assert len(xs & set(zt)) % 2 == 0


def test_small_bulk_commutes_and_counts():
    hx, hz, q = tile.build(T578, 4, 3, 3)
    assert len(q) == 2 * 6 * 6
    assert not ((hx.astype(np.int64) @ hz.T.astype(np.int64)) & 1).any()
    assert hx.sum(axis=1).max() == 8 and hz.sum(axis=1).max() == 8


def test_layout_is_one_qubit_per_site_within_the_bilayer_cap():
    hx, hz, q = tile.build(T578, 4, 5, 5)
    v = validate(tile.layout(q), check_supports(CSSCode(hx=hx, hz=hz, name="t", meta={})), layers=1)
    assert v["max_per_site"] == 1 and v["ok_sites"] and abs(v["radius"] - 37 ** 0.5) < 1e-9


@pytest.mark.skipif(not (BOARD / "578-18-20.json").exists(), reason="needs the upstream checkout in .worktrees/frontier")
def test_reproduces_board_578_18_20():
    hx, hz, q = tile.build(T578, 4, 14, 14)
    k = hx.shape[1] - gf2.rank(hx) - gf2.rank(hz)
    d = json.loads((BOARD / "578-18-20.json").read_text(encoding="utf-8"))

    def mat(rows):
        h = np.zeros((len(rows), d["n"]), dtype=np.uint8)
        for i, r in enumerate(rows):
            h[i, r] = 1
        return h
    board = CSSCode(hx=mat(d["checks"]["X"]), hz=mat(d["checks"]["Z"]), name="b", meta={})
    assert (hx.shape[1], k) == (578, 18)
    assert CSSCode(hx=hx, hz=hz, name="t", meta={}).signature() == board.signature()
