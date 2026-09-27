# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Every tracked code file in lab/ carries its SPDX header, and the licence files exist."""
from __future__ import annotations

import pathlib
import subprocess

LAB = pathlib.Path(__file__).resolve().parents[1]
HEADER = ("# SPDX-FileCopyrightText: 2026 Pete Salmond", "# SPDX-License-Identifier: Apache-2.0")


def _tracked_code() -> list[pathlib.Path]:
    out = subprocess.run(["git", "ls-files", "--", "."], cwd=LAB, capture_output=True, text=True, check=True).stdout
    return [LAB / f for f in out.split() if f.endswith((".py", ".sh"))]


def test_licence_files_exist():
    for name in ("LICENSE", "NOTICE", "CITATION.cff", "LICENSES/CC-BY-4.0.txt"):
        assert (LAB / name).is_file(), name
    assert "Apache License" in (LAB / "LICENSE").read_text(encoding="utf-8")
    assert "Attribution 4.0 International" in (LAB / "LICENSES/CC-BY-4.0.txt").read_text(encoding="utf-8")


def test_code_files_carry_spdx_header():
    files = _tracked_code()
    assert files, "git ls-files found no code under lab/"
    missing = []
    for f in files:
        head = f.read_text(encoding="utf-8").splitlines()[:4]
        if not all(h in head for h in HEADER):
            missing.append(str(f.relative_to(LAB)))
    assert not missing, f"missing SPDX header: {missing}"
