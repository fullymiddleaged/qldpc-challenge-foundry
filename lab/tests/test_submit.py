# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("submit", ROOT / "scripts/submit.py")
submit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(submit)

INFO = {"construction": "BB code", "layers": 2, "status": "submit"}


def test_refusal_blocks_hold_and_do_not_submit_unless_forced():
    for status in ("hold", "do_not_submit"):
        info = dict(INFO, status=status, status_reason="why")
        assert "refusing" in submit.refusal("x", info, force=False)
        assert submit.refusal("x", info, force=True) is None
    assert submit.refusal("x", INFO, force=False) is None
    assert submit.refusal("x", dict(INFO, status="optional"), force=False) is None


def test_qldpc_args_dry_run_with_layout():
    npz = Path("c.npz")
    cmd = submit.qldpc_args(INFO, npz, "@me", family="bivariate-bicycle", model="m", has_layout=True,
                            circuits=False, note=None, for_real=False)
    assert cmd[:2] == ["submit", "c.npz"]
    assert cmd[cmd.index("--coords") + 1] == "c.npz" and cmd[cmd.index("--layers") + 1] == "2"
    assert "--no-circuit" in cmd and cmd[-1] == "--dry-run" and "--open-pr" not in cmd


def test_qldpc_args_for_real_without_layout_or_model():
    cmd = submit.qldpc_args(INFO, Path("c.npz"), "@me", family="f", model="", has_layout=False,
                            circuits=True, note=Path("n.md"), for_real=True)
    assert "--coords" not in cmd and "--model" not in cmd and "--no-circuit" not in cmd
    assert cmd[cmd.index("--note-file") + 1] == "n.md"
    assert cmd[-1] == "--open-pr" and "--dry-run" not in cmd


def test_qldpc_args_json_only_on_dry_run():
    kw = dict(family="f", model="m", has_layout=False, circuits=False, note=None)
    assert submit.qldpc_args(INFO, Path("c.npz"), "@me", **kw, for_real=False, json_doc=True)[-1] == "--json"
    assert "--json" not in submit.qldpc_args(INFO, Path("c.npz"), "@me", **kw, for_real=True, json_doc=True)


def test_doc_from_dry_run_parses_the_json_after_the_marker():
    out = "  OK  verified. score {x}\n\n--dry-run: would write codes/1-1-1.json\n{\n \"n\": 1,\n \"k\": 1\n}\n"
    assert submit.doc_from_dry_run(out) == {"n": 1, "k": 1}
    try:
        submit.doc_from_dry_run("verification FAILED")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")


def test_gate_verdict():
    assert submit.gate_verdict({"passed": True, "labels": []}) is None
    assert "refuted" in submit.gate_verdict({"passed": False, "labels": ["refuted (over-claimed distance)"]})
    assert submit.gate_verdict({"passed": False}).endswith("no labels")


GENOME = {"l": 6, "m": 6, "A": [[3, 0], [0, 1], [0, 2]], "B": [[0, 3], [1, 0], [2, 0]]}
CAND = dict(INFO, params="[[72,12,6]]", genome=GENOME)


def _write_proof(tmp_path, **over):
    proof = {"result": "minimum distance proven exact", "d": 6, "genome": GENOME, **over}
    (tmp_path / "distance_proof_72_12_6.json").write_text(submit.json.dumps(proof), encoding="utf-8")


def test_params_parses_the_manifest_triple():
    assert submit.params({"params": "[[216,8,16]]"}) == (216, 8, 16)


def test_proof_refusal_accepts_a_matching_finished_proof(tmp_path):
    # same genome with the terms listed in another order is the same code
    _write_proof(tmp_path, genome=dict(GENOME, A=[[0, 2], [3, 0], [0, 1]]))
    assert submit.proof_refusal("c", CAND, tmp_path) is None


def test_proof_refusal_without_a_proof_file(tmp_path):
    assert "no distance proof" in submit.proof_refusal("c", CAND, tmp_path)


def test_proof_refusal_on_an_unfinished_proof(tmp_path):
    _write_proof(tmp_path, result="NOT YET PROVEN. Upper bound 6")
    assert "not a finished proof" in submit.proof_refusal("c", CAND, tmp_path)


def test_proof_refusal_on_a_different_distance(tmp_path):
    _write_proof(tmp_path, d=5)
    assert "proves d=5, the manifest claims d=6" in submit.proof_refusal("c", CAND, tmp_path)


def test_proof_refusal_on_a_different_genome(tmp_path):
    _write_proof(tmp_path, genome=dict(GENOME, B=[[0, 3], [1, 0], [3, 0]]))
    assert "different genome" in submit.proof_refusal("c", CAND, tmp_path)


def test_claim_mismatch():
    assert submit.claim_mismatch({"distance": {"d": 16}}, 16) is None
    assert submit.claim_mismatch({"distance": {"d": 18}}, 16) == "qldpc would submit d=18, but the proven distance is 16"


def test_real_candidates_marked_submit_have_matching_proofs():
    manifest = submit.json.loads((ROOT / "results/candidates/manifest.json").read_text(encoding="utf-8"))
    for name, info in manifest.items():
        if info.get("status") == "submit":
            assert submit.proof_refusal(name, info, submit.PROOFS) is None, name


def test_qldpc_args_passes_provenance_notes_only_when_set():
    kw = dict(family="f", model="m", has_layout=False, circuits=False, note=None, for_real=False)
    cmd = submit.qldpc_args(dict(INFO, provenance_notes="Checked, not equivalent."), Path("c.npz"), "@me", **kw)
    assert cmd[cmd.index("--notes") + 1] == "Checked, not equivalent."
    assert "--notes" not in submit.qldpc_args(INFO, Path("c.npz"), "@me", **kw)


def test_free_slug_takes_the_next_suffix(tmp_path):
    assert submit.free_slug(tmp_path, 144, 12, 12) == "144-12-12"
    (tmp_path / "144-12-12.json").write_text("{}")
    assert submit.free_slug(tmp_path, 144, 12, 12) == "144-12-12-b"
    (tmp_path / "144-12-12-b.json").write_text("{}")
    assert submit.free_slug(tmp_path, 144, 12, 12) == "144-12-12-c"


def test_proof_refusal_handles_tile_genomes(tmp_path):
    tile = {"tile": [["h", 0, 0], ["v", 1, 2]], "B": 3, "l": 11, "m": 11}
    info = {"params": "[[338,8,16]]", "genome": tile}
    proof = {"result": "minimum distance proven exact", "d": 16, "genome": {**tile, "tile": [["v", 1, 2], ["h", 0, 0]]}}
    (tmp_path / "distance_proof_338_8_16.json").write_text(json.dumps(proof))
    assert submit.proof_refusal("t", info, tmp_path) is None                       # edge order does not matter
    proof["genome"] = {**tile, "l": 12}
    (tmp_path / "distance_proof_338_8_16.json").write_text(json.dumps(proof))
    assert "different genome" in submit.proof_refusal("t", info, tmp_path)
