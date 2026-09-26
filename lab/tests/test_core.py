"""Core regression tests. Run with `pytest -q`, or without pytest via
`python tests/run_all.py`. Everything here is pure numpy/scipy (no stim or
ldpc needed) and finishes in about a minute."""
from __future__ import annotations

import itertools
import json
import tempfile
from pathlib import Path

import numpy as np

from qec_search import gf2
from qec_search.bbcode import BBGenome, build_bb
from qec_search.known_codes import KNOWN, KNOWN_SIGS

ROOT = Path(__file__).resolve().parents[1]
X = lambda i: (i, 0)
Y = lambda j: (0, j)


def test_known_codes_reproduce_published_parameters():
    for name, g, (n, k, d) in KNOWN:
        c = build_bb(g)
        assert c.check_css(), name
        assert (c.n, c.k) == (n, k), name
        # upper bound must reach the published distance and never undershoot it
        assert c.distance_upper_bound(trials=400, stop_at=d) == d, name


def test_signature_identifies_relabelled_gross_code():
    gross = build_bb(KNOWN[3][1])
    variant = build_bb(BBGenome(12, 6, (X(9), Y(1), Y(2)), (Y(3), X(1), X(2))))
    assert gross.signature() == variant.signature()
    assert variant.signature() in KNOWN_SIGS
    other = build_bb(BBGenome(12, 6, (Y(1), X(3), X(10)), (Y(3), X(1), X(2))))  # [[144,8,12]]
    assert other.signature() != gross.signature()


def test_split_code_is_detected():
    two_gross = build_bb(BBGenome(24, 6, (Y(3), X(2), X(4)), (Y(1), Y(2), X(6))))  # [[288,24,12]]
    assert two_gross.k == 24 and two_gross.components() == 2
    assert build_bb(KNOWN[3][1]).components() == 1


def test_exact_distance_milp_small_codes():
    from qec_search.exact_distance import certify_bb
    for idx, d in ((0, 6), (1, 10)):  # [[72,12,6]], [[90,8,10]]
        g = KNOWN[idx][1]
        r = certify_bb(build_bb(g), g.l, g.m, upper=None, time_limit=120)
        assert r["exact"] and r["best_weight"] == d


def _independent_candidate_check(npz_path):
    with np.load(npz_path) as z:
        hx, hz, coords = z["hx"], z["hz"], z["coords"]
    n = hx.shape[1]
    assert not ((hx.astype(int) @ hz.T.astype(int)) % 2).any()
    k = n - gf2.rank(hx) - gf2.rank(hz)
    w = int(max(hx.sum(1).max(), hz.sum(1).max()))
    radius = 0.0
    for row in np.vstack([hx, hz]):
        s = np.nonzero(row)[0]
        P = coords[s]
        radius = max(radius, max(float(np.hypot(*(P[a] - P[b]))) for a, b in itertools.combinations(range(len(s)), 2)))
    sites, counts = np.unique(coords, axis=0, return_counts=True)
    gaps = np.sqrt(((sites[:, None] - sites[None]) ** 2).sum(-1)) + np.eye(len(sites)) * 99
    return n, k, w, radius, int(counts.max()), float(gaps.min())


def test_record_candidates_meet_bilayer_rules():
    manifest = json.load(open(ROOT / "results/candidates/manifest.json"))
    for name, info in manifest.items():
        n, k, w, radius, per_site, gap = _independent_candidate_check(ROOT / f"results/candidates/{name}.npz")
        pn, pk, _ = (int(x) for x in info["params"].strip("[]").split(","))
        assert (n, k) == (pn, pk), name
        assert w <= 6 and radius <= 7.0 and per_site <= 2 and gap >= 1.0, name


def test_layout_finder_on_small_code():
    from qec_search.layout import find_layout, validate, check_supports
    g = KNOWN[3][1]  # gross code fits quickly
    c = build_bb(g)
    coords, dmax = find_layout(c, g, steps=300_000, restarts=2)
    v = validate(coords, check_supports(c))
    assert v["ok_sites"] and v["radius"] <= 7.0


def test_numpy_decoder_corrects_errors():
    from qec_search.evaluate import code_capacity_ler
    c = build_bb(KNOWN[0][1])
    r = code_capacity_ler(c, 0.01, shots=1000, max_fails=1000, force_numpy=True)
    assert r["ler"] < 0.05


def test_export_roundtrip():
    from qec_search import export
    with tempfile.TemporaryDirectory() as tmp:
        run = Path(tmp) / "run"
        run.mkdir()
        recs = []
        for g, d in ((BBGenome(12, 6, (Y(1), X(3), X(10)), (Y(3), X(1), X(2))), 12),
                     (KNOWN[3][1], 12)):  # the second is IBM's and must be skipped
            c = build_bb(g)
            recs.append({"key": g.canonical_key(), "genome": g.to_json(), "pretty": g.pretty(), "n": c.n,
                         "k": c.k, "d_ub": d, "kd2n": c.k * d * d / c.n})
        (run / "results.jsonl").write_text("\n".join(json.dumps(r) for r in recs) + "\n")
        dest = Path(tmp) / "out"
        export.main([str(run), "--dest", str(dest)])
        m = json.load(open(dest / "manifest.json"))
        assert list(m) == ["bb_144_8_12_12x6"]
        with np.load(dest / "bb_144_8_12_12x6.npz") as z:  # close it, or Windows can't delete tmp
            assert z["hx"].shape == (72, 144)
