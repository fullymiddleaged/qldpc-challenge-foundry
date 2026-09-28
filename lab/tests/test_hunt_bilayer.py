# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""scripts/hunt_bilayer.py: frontier parsing, Pareto dominance, and covering-graph lifts of BB codes.

The cover test uses the property from arXiv:2511.13560 that an h-fold cover of a BB code keeps at least its k, and
that the trivial lift (every t = 0) is among the covers. It runs on the gross code."""
from __future__ import annotations

import importlib.util
import pathlib

from qec_search.bbcode import BBGenome, build_bb

LAB = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("hunt_bilayer", LAB / "scripts" / "hunt_bilayer.py")
hunt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hunt)

GROSS = BBGenome(12, 6, ((3, 0), (0, 1), (0, 2)), ((0, 3), (1, 0), (2, 0)))
TARGETS = """2D-local bilayer / weight <= 6
  387 codes, 253 nondominated, best kd2/n 19.20
    [[288,12,16]] w=6 kd2/n=10.67
    [[300,8,20]] w=6 kd2/n=10.67
"""


def test_parse_frontier_reads_every_entry():
    assert hunt.parse_frontier(TARGETS) == [(288, 12, 16, 6), (300, 8, 20, 6)]


def test_dominance_ties_and_each_axis():
    f = hunt.parse_frontier(TARGETS)
    assert hunt.dominated((288, 12, 16, 6), f)          # a tie is no record
    assert hunt.dominated((300, 12, 14, 6), f)          # worse n and d
    assert not hunt.dominated((288, 8, 20, 6), f)       # beats [[300,8,20]] on n, loses to nothing
    assert not hunt.dominated((288, 12, 17, 6), f)      # better d
    assert not hunt.dominated((286, 12, 16, 6), f)      # better n
    assert hunt.dominated((288, 12, 16, 8), f)          # heavier checks lose the tie


def test_double_covers_of_the_gross_code_keep_k():
    base_k = build_bb(GROSS).k
    cs = hunt.covers(GROSS, 2)
    assert {(c.l, c.m) for c in cs} == {(24, 6), (12, 12)}
    assert GROSS.canonical_key() != cs[0].canonical_key()
    trivial = BBGenome(24, 6, GROSS.A, GROSS.B).canonical_key()
    assert trivial in {c.canonical_key() for c in cs}
    assert all(build_bb(c).k >= base_k for c in cs)


def test_screen_keeps_what_nothing_beats_and_drops_ties():
    kept = hunt.screen_one(GROSS.to_json(), [], kmin=4, trials=20)
    assert kept is not None and (kept["n"], kept["k"], kept["w"]) == (144, 12, 6) and kept["d_ub"] >= 12
    assert hunt.screen_one(GROSS.to_json(), [(144, 12, 12, 6)], kmin=4, trials=20) is None
    assert hunt.screen_one(GROSS.to_json(), [], kmin=13, trials=20) is None
