# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""scripts/audit_board_distances.py: reading DistQLDPC logs and turning per-side results into a verdict."""
from __future__ import annotations

import importlib.util
import pathlib

LAB = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("audit", LAB / "scripts" / "audit_board_distances.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

FINISHED = "c trying d: 8\nc d_lb: 8\nc d_ub: 14\nc d_ub: 12\nc d_lb: 12\nc d  : 12\no 12\nwall 18.9 s\n"
TIMED_OUT = "c d_lb: 8\nc d_ub: 16\nc d_lb: 11\nc d_ub: 15\nwall 900.1 s\n"


def side(exact=None, lb=None, ub=None):
    return {"exact": exact, "lb": lb, "ub": ub}


def test_parse_finished_and_timed_out_logs():
    assert audit.parse_log(FINISHED) == {"exact": 12, "lb": 12, "ub": 12}
    assert audit.parse_log(TIMED_OUT) == {"exact": None, "lb": 11, "ub": 15}
    assert audit.parse_log("") == {"exact": None, "lb": None, "ub": None}


def test_verdicts():
    assert audit.verdict(12, {"X": side(12), "Z": side(13)}) == "holds"
    assert audit.verdict(12, {"X": side(11), "Z": side(12)}) == "LOWER"
    assert audit.verdict(12, {"X": side(None, 9, 11), "Z": side(12)}) == "LOWER"    # an unfinished side can refute
    assert audit.verdict(12, {"X": side(None, 10, 12), "Z": side(12)}) == "open"
    assert audit.verdict(12, {"X": side(None, 12, 12), "Z": side(None, 12, 13)}) == "holds"
    assert audit.verdict(12, {"X": side(13), "Z": side(14)}) == "ERROR"
