# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Tanner-graph automorphisms, orbital branching and the symmetry-reduced SAT certifier.

The end-to-end cases use small board codes whose distance upstream certified exactly. Each must come out PROVEN at
its d and REFUTED, with a valid witness, at d + 1. [[37,1,7]] is the regression for upstream's incomplete
sat_certify._logicals basis, which "proves" d = 8 on that code."""
from __future__ import annotations

import importlib.util
import pathlib
import tempfile

import numpy as np
import pytest

from qec_search import certify_sym
from qec_search.automorphisms import TannerGraph, find_duality, maps_checks, qubit_orbits

LAB = pathlib.Path(__file__).resolve().parents[1]
CODES = LAB.parent / "codes"
GROSS = LAB / "results" / "candidates" / "ibm_gross_144_12_12.npz"
HAS_SAT = all(importlib.util.find_spec(m) for m in ("pycryptosat", "pysat"))


def _gross():
    z = np.load(GROSS)
    return z["hx"].astype(np.int8), z["hz"].astype(np.int8)


def test_gross_code_is_transitive_and_self_dual():
    hx, hz = _gross()
    orbits, gens = qubit_orbits(hx, hz)
    assert [len(o) for o in orbits] == [144]
    g = TannerGraph(hx, hz)
    assert gens and all(maps_checks(g, g, p) for p in gens)
    assert find_duality(hx, hz) is not None


def test_bogus_permutation_is_rejected():
    hx, hz = _gross()
    g = TannerGraph(hx, hz)
    p = np.arange(144)
    p[[0, 1]] = p[[1, 0]]
    assert not maps_checks(g, g, p)


def test_stabiliser_orbits_respect_among():
    hx, hz = _gross()
    free = set(range(1, 72))
    orbits, gens = qubit_orbits(hx, hz, fixed=(0,), among=free)
    assert sorted(q for o in orbits for q in o) == sorted(free)
    assert all(set(p[sorted(free)].tolist()) == free and p[0] == 0 for p in gens)


def test_breaking_the_symmetry_shrinks_the_group():
    hx, hz = _gross()
    hx = hx.copy()
    hx[0, np.nonzero(hx[0])[0][0]] ^= 1                  # perturb one X check
    orbits, _ = qubit_orbits(hx, hz)
    assert len(orbits) > 1


def test_orbital_cubes_are_disjoint_on_their_heads():
    hx, hz, _ = certify_sym.load_code(CODES / "56-6-7.json")
    cubes = certify_sym.orbital_cubes(hx, hz, depth=2)
    heads = [tuple(sorted(x for x in c if x > 0)) for c in cubes]
    assert len(set(heads)) == len(heads) and len(cubes) > 1


@pytest.mark.skipif(not HAS_SAT, reason="needs pycryptosat and python-sat (pip install -e .[cert])")
def test_pairing_set_is_complete_and_upstream_now_agrees():
    hx, hz, _ = certify_sym.load_code(CODES / "37-1-7.json")
    gf2 = certify_sym._upstream().gf2
    t = certify_sym.pairing_set(hx, hz)
    assert len(t) == 1 and gf2.rank(np.vstack([hz, t])) == gf2.rank(hz) + 1
    # upstream's _logicals returned a stabiliser here (no logical class) until our issue #2273, fixed by #2343
    upstream = certify_sym._upstream()._logicals(hz, hx)
    assert gf2.rank(np.vstack([hz, upstream])) == gf2.rank(hz) + 1


@pytest.mark.skipif(not HAS_SAT, reason="needs pycryptosat and python-sat (pip install -e .[cert])")
@pytest.mark.parametrize("name", ["37-1-7", "30-4-6", "54-6-9", "56-6-7"])
def test_certifier_proves_d_and_refutes_d_plus_one(name, monkeypatch):
    monkeypatch.setattr(certify_sym, "LOGS", pathlib.Path(tempfile.mkdtemp()))
    path = CODES / f"{name}.json"
    hx, hz, d = certify_sym.load_code(path)
    for depth in (0, 2):
        st = certify_sym.certify([path], d, depth, True, 120, 2, tag=f"_t{depth}", quiet=True)[path]
        assert certify_sym.summarize(st)["verdict"] == "PROVEN"
        st = certify_sym.certify([path], d + 1, depth, True, 120, 2, tag=f"_t{depth}", quiet=True)[path]
        assert certify_sym.summarize(st)["verdict"] == "REFUTED"
        wit = next(c["witness"] for c in st["cubes"].values() if c["status"] == "SAT")
        side = next(k.split("|")[0] for k, c in st["cubes"].items() if c["status"] == "SAT")
        h_same, h_opp = (hx, hz) if side == "X" else (hz, hx)
        assert len(wit) <= d and certify_sym.is_logical(h_same, h_opp, wit)


def _min_weight_logicals(h_same, h_opp, d):
    """Every logical in ker(h_opp) outside rowspace(h_same) of weight exactly d, by DFS with incremental syndrome."""
    n = h_opp.shape[1]
    t = certify_sym.pairing_set(h_same, h_opp)
    col = [int("".join(map(str, h_opp[::-1, q])), 2) for q in range(n)]
    tcol = [int("".join(map(str, t[::-1, q])), 2) for q in range(n)]
    out = []

    def dfs(start, depth, syn, anti, chosen):
        if depth == d:
            if syn == 0 and anti:
                out.append(frozenset(chosen))
            return
        for q in range(start, n - (d - depth) + 1):
            chosen.append(q)
            dfs(q + 1, depth + 1, syn ^ col[q], anti ^ tcol[q], chosen)
            chosen.pop()

    dfs(0, 0, 0, 0, [])
    return set(out)


def _orbits_missed(name, cubes_fn, depth):
    """How many Aut-orbits of minimum-weight X logicals have no member inside any root cube (0 = sound)."""
    hx, hz, d = certify_sym.load_code(CODES / f"{name}.json")
    logicals = _min_weight_logicals(hx, hz, d)
    _, gens = qubit_orbits(hx, hz)
    cubes = cubes_fn(hx, hz, depth)
    seen, missed = set(), 0
    for s in logicals:
        if s in seen:
            continue
        orb, stack = {s}, [s]
        while stack:
            x = stack.pop()
            for g in gens:
                y = frozenset(int(g[q]) for q in x)
                assert y in logicals, "an automorphism moved a minimum-weight logical outside the set"
                if y not in orb:
                    orb.add(y)
                    stack.append(y)
        seen |= orb
        if not any(all((lit > 0) == (abs(lit) - 1 in s) for lit in c) for s in orb for c in cubes):
            missed += 1
    return missed


@pytest.mark.skipif(not HAS_SAT, reason="needs pycryptosat and python-sat (pip install -e .[cert])")
@pytest.mark.parametrize("name", ["45-9-3", "37-1-7"])
def test_orbital_branching_keeps_every_orbit_of_lightest_logicals(name):
    for depth in (1, 2, 3):
        assert _orbits_missed(name, certify_sym.orbital_cubes, depth) == 0


@pytest.mark.skipif(not HAS_SAT, reason="needs pycryptosat and python-sat (pip install -e .[cert])")
def test_coverage_check_catches_a_dropped_cube():
    def broken(hx, hz, depth):
        return certify_sym.orbital_cubes(hx, hz, depth)[1:]
    assert _orbits_missed("45-9-3", broken, 2) > 0
