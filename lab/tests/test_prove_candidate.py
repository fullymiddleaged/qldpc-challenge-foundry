# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
from __future__ import annotations

import importlib.util
import pathlib

_spec = importlib.util.spec_from_file_location("prove_candidate", pathlib.Path(__file__).parents[1] / "scripts" / "prove_candidate.py")
pc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pc)

INFO = {"params": "[[360,8,17]]", "genome": {"tile": [], "B": 3, "l": 12, "m": 12}}


def test_both_sides_exact_at_the_claimed_d_is_a_proof():
    sides = {"X": {"exact": 17, "lb": 17, "ub": 17}, "Z": {"exact": 18, "lb": 18, "ub": 18}}
    rec = pc.proof_record("c", INFO, sides, 60, "2026-10-16")
    assert rec["result"] == "minimum distance proven exact" and rec["d"] == 17 and rec["code"] == "[[360,8,17]]"


def test_exact_but_lighter_than_claimed_is_flagged():
    sides = {"X": {"exact": 16, "lb": 16, "ub": 16}, "Z": {"exact": 17, "lb": 17, "ub": 17}}
    assert pc.proof_record("c", INFO, sides, 60, "d")["result"] == "proven d = 16, manifest claims 17"


def test_an_unfinished_side_is_open_with_its_bound():
    sides = {"X": {"exact": 17, "lb": 17, "ub": 17}, "Z": {"exact": None, "lb": 15, "ub": 17}}
    rec = pc.proof_record("c", INFO, sides, 60, "d")
    assert rec["result"].startswith("open: proven d >= 15") and "proven exact" not in rec["result"]
    assert pc.proof_record("c", INFO, {"X": {"exact": 17}}, 60, "d")["result"].startswith("open")
