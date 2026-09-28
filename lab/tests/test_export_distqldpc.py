# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""scripts/export_distqldpc.py poses one CSS side for DistQLDPC: the z variables must satisfy exactly the X-side
conditions of certify_sym (in ker Hz, anticommuting with a complete logical basis), and a cube's zeros become unit
rows. Checked on the gross code, without the solver itself."""
from __future__ import annotations

import importlib.util
import pathlib
import tempfile

import numpy as np
import pytest

from qec_search import certify_sym

LAB = pathlib.Path(__file__).resolve().parents[1]
GROSS = LAB / "results" / "candidates" / "ibm_gross_144_12_12.npz"
HAS_GF2 = (LAB.parent / "verify" / "sat_certify.py").exists()

spec = importlib.util.spec_from_file_location("export_distqldpc", LAB / "scripts" / "export_distqldpc.py")
export = importlib.util.module_from_spec(spec)
spec.loader.exec_module(export)


def _gross():
    z = np.load(GROSS)
    return z["hx"].astype(np.int8), z["hz"].astype(np.int8)


def _anticommutes(gz: np.ndarray, v: np.ndarray) -> bool:
    return bool((gz.astype(np.int64) @ v % 2).any())


@pytest.mark.skipif(not HAS_GF2, reason="needs upstream verify/sat_certify.py")
def test_x_side_accepts_logicals_and_rejects_stabilisers():
    hx, hz = _gross()
    m = export.one_sided(hx, hz, "X")
    assert m["Hz"].shape == (144, 144) and not m["Gx"].any()
    for x in certify_sym.pairing_set(hz, hx):                    # X-type logicals, one per class
        assert not (m["Hx"].astype(np.int64) @ x % 2).any()
        assert _anticommutes(m["Gz"], x)
    for s in hx:                                                 # X stabilisers are trivial
        assert not _anticommutes(m["Gz"], s)


@pytest.mark.skipif(not HAS_GF2, reason="needs upstream verify/sat_certify.py")
def test_cube_zeros_become_unit_rows_and_files_parse():
    hx, hz = _gross()
    m = export.one_sided(hx, hz, "X", [3, 7])
    extra = m["Hx"][len(hz):]
    assert extra.shape == (2, 144) and list(np.nonzero(extra)[1]) == [3, 7] and extra.sum() == 2
    with tempfile.TemporaryDirectory() as d:
        prefix = pathlib.Path(d) / "gross_X"
        export.write(m, prefix)
        for name, mat in m.items():
            text = pathlib.Path(f"{prefix}_{name}.txt").read_bytes()
            assert b"\r" not in text                             # read by a Linux build
            rows = [[int(c) for c in line.split()] for line in text.decode().splitlines()]
            assert np.array_equal(np.array(rows, dtype=np.int8), mat)
