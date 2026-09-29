# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Exact CSS distances with DistQLDPC (the MaxCDCL solver of arXiv:2606.12445), the one path every lab tool uses.

Every speed-up lives here so the hunt, the board audit and anything later get it together:
  * one side at a time (Hz := I, Gx := a zero row), about 4.6x faster than the solver's two-sided form;
  * X/Z duality (automorphisms.find_duality): one side suffices when it exists;
  * orbital branching (certify_sym.orbital_cubes): each side splits into cubes, one per qubit orbit of the Tanner
    graph's automorphism group; on the 180 coset codes this took a distance proof from over 75 cpu-h to minutes.
The minimum over a side's cubes is that side's distance, because the cubes cover every logical up to automorphism.

The solver is built in WSL (~/src/DistQLDPC, with scripts/distqldpc-one-z.patch) and run through
scripts/run_distqldpc.sh; logs land in results/distqldpc/logs/<stem>.log.
"""
from __future__ import annotations

import os
import pathlib
import re
import subprocess
from dataclasses import dataclass

import numpy as np

LAB = pathlib.Path(__file__).resolve().parents[1]
DQ = LAB / "results" / "distqldpc"
WSL_LAB = "/mnt/c/vscode/qldpc-challenge-foundry/lab"


@dataclass(frozen=True)
class Job:
    code: str        # caller's name for the code
    side: str        # "X": an X-type logical, x in ker(Hz)
    cube: tuple      # certify_sym cube: 1-based literals, positive = forced into the support
    stem: str        # file stem under results/distqldpc

    @property
    def flag(self) -> str:
        ones = [v - 1 for v in self.cube if v > 0]
        return f"-one-z={','.join(map(str, ones))}" if ones else ""


def one_sided(hx: np.ndarray, hz: np.ndarray, side: str, zeros=()) -> dict[str, np.ndarray]:
    """DistQLDPC input posing one CSS side; zeros (0-based qubits) become unit rows of Hx."""
    from .certify_sym import pairing_set
    h_same, h_opp = (hx, hz) if side == "X" else (hz, hx)
    n = hx.shape[1]
    fixed = np.zeros((len(zeros), n), dtype=np.int8)
    fixed[np.arange(len(zeros)), list(zeros)] = 1
    return {"Hx": np.vstack([h_opp, fixed]), "Hz": np.eye(n, dtype=np.int8), "Gx": np.zeros((1, n), dtype=np.int8),
            "Gz": pairing_set(h_same, h_opp)}


def write(mats: dict[str, np.ndarray], prefix: pathlib.Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for name, m in mats.items():
        lines = [" ".join(map(str, row)) for row in np.asarray(m, dtype=np.int8)]
        pathlib.Path(f"{prefix}_{name}.txt").write_text("\n".join(lines) + "\n", newline="\n")


def plan(name: str, hx: np.ndarray, hz: np.ndarray, depth: int = 1, max_cubes: int = 16) -> list[Job]:
    """Sides (one if a duality exists) times orbital cubes, for one code. A code with few automorphisms has many
    orbits, and one capped solver run per orbit is slower than one unsplit run (a [[216,4]] hunt code split into 169
    cubes on 29 Sep), so past max_cubes the side is solved whole."""
    from .automorphisms import find_duality
    from .certify_sym import orbital_cubes
    sides = ["X"] if find_duality(hx, hz) is not None else ["X", "Z"]
    cubes = orbital_cubes(hx, hz, depth)
    if len(cubes) > max_cubes:
        cubes = [[]]
    return [Job(name, s, tuple(c), f"{name}_{s}_c{i}") for s in sides for i, c in enumerate(cubes)]


def parse_log(text: str) -> dict:
    """Optimum if the solver finished, else the last bounds it reported."""
    m = re.search(r"^o (\d+)", text, re.M)
    lbs = [int(x) for x in re.findall(r"^c d_lb: (\d+)", text, re.M)]
    ubs = [int(x) for x in re.findall(r"^c d_ub: (\d+)", text, re.M)]
    return {"exact": int(m.group(1)) if m else None, "lb": max(lbs, default=None), "ub": min(ubs, default=None)}


def combine(cubes: list[dict]) -> dict:
    """One side from its cubes: exact when every cube finished; else the proven bound (min over cubes of what each
    proved) and the best logical found in any cube."""
    if cubes and all(c["exact"] is not None for c in cubes):
        e = min(c["exact"] for c in cubes)
        return {"exact": e, "lb": e, "ub": e}
    lb = min((c["exact"] if c["exact"] is not None else (c["lb"] or 0)) for c in cubes)
    ubs = [c["exact"] if c["exact"] is not None else c["ub"] for c in cubes]
    ubs = [u for u in ubs if u is not None]
    return {"exact": None, "lb": lb, "ub": min(ubs) if ubs else None}


def solve_many(codes: dict[str, tuple[np.ndarray, np.ndarray]], secs: int, jobs: int = os.cpu_count(),
               depth: int = 1) -> dict[str, dict]:
    """{name: (hx, hz)} -> {name: {"sides": {side: {exact, lb, ub}}, "cubes": {stem: {...}}}}, all in one WSL batch.
    A side is exact only when all its cubes finished within secs."""
    all_jobs = [j for name, (hx, hz) in codes.items() for j in plan(name, hx, hz, depth)]
    for j in all_jobs:
        hx, hz = codes[j.code]
        write(one_sided(hx, hz, j.side, [-v - 1 for v in j.cube if v < 0]), DQ / j.stem)
    args = [f"results/distqldpc/{j.stem}" + (f":{j.flag}" if j.flag else "") for j in all_jobs]
    env = dict(os.environ, JOBS=str(jobs), WSLENV="JOBS", MSYS_NO_PATHCONV="1")
    subprocess.run(["wsl.exe", "-d", "Ubuntu", "--", "bash", f"{WSL_LAB}/scripts/run_distqldpc.sh", str(secs), *args],
                   env=env, capture_output=True)
    out: dict[str, dict] = {name: {"sides": {}, "cubes": {}} for name in codes}
    per_side: dict[tuple[str, str], list[dict]] = {}
    for j in all_jobs:
        log = DQ / "logs" / f"{j.stem}.log"
        r = parse_log(log.read_text() if log.exists() else "")
        out[j.code]["cubes"][j.stem] = r
        per_side.setdefault((j.code, j.side), []).append(r)
    for (name, side), rs in per_side.items():
        out[name]["sides"][side] = combine(rs)
    return out
