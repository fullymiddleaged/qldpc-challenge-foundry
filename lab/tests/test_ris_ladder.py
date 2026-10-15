# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
from __future__ import annotations

import importlib.util
import json
import pathlib

_spec = importlib.util.spec_from_file_location("ris_ladder", pathlib.Path(__file__).parents[1] / "scripts" / "ris_ladder.py")
ladder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ladder)


def test_done_rungs_reads_the_log_so_a_rerun_resumes(tmp_path):
    log = tmp_path / "ladder.jsonl"
    assert ladder.done_rungs(log) == set()
    log.write_text("\n".join(json.dumps({"npz": "cand/a.npz", "trials": t, "d_ub": 17}) for t in (20000, 200000)) + "\n\n")
    assert ladder.done_rungs(log) == {("cand/a.npz", 20000), ("cand/a.npz", 200000)}


def test_rung_seeds_are_stable_and_distinct_per_code_and_rung():
    assert ladder.rung_seed("cand/a.npz", 20000) == ladder.rung_seed("cand/a.npz", 20000)
    seeds = {ladder.rung_seed(f"cand/{c}.npz", t) for c in "ab" for t in (20000, 200000, 1000000)}
    assert len(seeds) == 6 and all(0 <= s < 2**31 for s in seeds)
