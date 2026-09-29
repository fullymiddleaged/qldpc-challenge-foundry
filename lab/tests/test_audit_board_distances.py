# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""scripts/audit_board_distances.py and qec_search.distqldpc: reading solver logs, combining orbital cubes into a
side, and turning per-side results into a verdict."""
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


def test_combine_cubes_into_a_side():
    from qec_search.distqldpc import combine
    assert combine([side(12), side(36)]) == {"exact": 12, "lb": 12, "ub": 12}
    assert combine([side(None, 9, 14), side(20)]) == {"exact": None, "lb": 9, "ub": 14}    # unfinished cube bounds it
    assert combine([side(None, 11, None), side(13)]) == {"exact": None, "lb": 11, "ub": 13}


def test_plan_uses_duality_and_orbits_on_the_gross_code():
    import numpy as np
    from qec_search.distqldpc import plan
    z = np.load(LAB / "results" / "candidates" / "ibm_gross_144_12_12.npz")
    jobs = plan("gross", z["hx"].astype(np.int8), z["hz"].astype(np.int8))
    assert {j.side for j in jobs} == {"X"}                 # d_X = d_Z, so one side
    assert len(jobs) == 1 and jobs[0].flag == "-one-z=0"   # one qubit orbit: one cube, qubit 0 forced in
