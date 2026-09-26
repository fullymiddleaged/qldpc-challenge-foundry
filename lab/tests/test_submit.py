from __future__ import annotations

import importlib.util
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
