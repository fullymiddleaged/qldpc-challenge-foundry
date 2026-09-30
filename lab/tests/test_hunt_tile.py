# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""scripts/hunt_tile.py: the tile radius filter, the tile enumeration, and the screen's frontier test."""
from __future__ import annotations

import importlib.util
import math
import pathlib

LAB = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("hunt_tile", LAB / "scripts" / "hunt_tile.py")
ht = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ht)

T578 = [("h", 0, 0), ("h", 0, 3), ("h", 2, 2), ("h", 3, 0), ("v", 0, 1), ("v", 1, 1), ("v", 2, 0), ("v", 3, 3)]


def test_radius_matches_the_board_entry():
    assert math.isclose(ht.tile_radius(T578, 4), 37 ** 0.5)       # codes/578-18-20: interaction radius sqrt(37)


def test_exhaustive_enumeration_size():
    assert sum(1 for _ in ht.all_tiles(2, 4, 0, 0)) == math.comb(8, 4)
    assert sum(1 for _ in ht.all_tiles(4, 6, 25, 0)) == 25             # random above B = 3


def test_screen_respects_the_cap_and_the_frontier():
    small = [("h", 0, 0), ("h", 0, 1), ("h", 1, 0), ("v", 0, 0)]      # connected, k = 2 at 6 x 6
    assert ht.screen_one(T578, 4, 8, [], 2, 200, 20, 4.0) == []       # radius 6.08 > 4: never built
    kept = ht.screen_one(small, 2, 4, [], 1, 120, 20, 4.0)
    assert kept and all(r["n"] <= 120 and r["radius"] <= 4.0 for r in kept)
    beaten = [(1, 10 ** 3, 10 ** 3, 4)]                                # an entry that beats everything
    assert ht.screen_one(small, 2, 4, beaten, 1, 120, 20, 4.0) == []
