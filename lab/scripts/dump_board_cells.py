# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Load the board once with upstream's own site/build.py and write each 2D-local cell's nondominated entries.

`qldpc targets` re-verifies every code on each call (about 45 min here); this does it once and writes
results/hunt/board_entries.json plus results/hunt/frontier_<locality>_w<weight>.txt in the '[[n,k,d]] w=W' form
scripts/hunt_bilayer.py reads. Run it from an upstream checkout (e.g. .worktrees/frontier at upstream/main):

    PYTHONUTF8=1 uv run python ../../lab/scripts/dump_board_cells.py
"""
from __future__ import annotations

import json
import pathlib
import sys

LAB = pathlib.Path(__file__).resolve().parents[1]
OUT = LAB / "results" / "hunt"


def main() -> None:
    sys.path.insert(0, str(pathlib.Path.cwd() / "site"))
    from build import cell_key, cells, load_entries, pareto
    entries = load_entries()
    OUT.mkdir(parents=True, exist_ok=True)
    keep = ("slug", "n", "k", "d", "w", "code_type", "locality_class", "weight_class")
    (OUT / "board_entries.json").write_text(json.dumps([{k: e.get(k) for k in keep} for e in entries], indent=0))
    by_cell = {}
    for e in entries:
        for cell in cells(e):
            by_cell.setdefault(cell_key(e, cell), []).append(e)
    names = {"local-2d-bilayer": "bilayer", "local-2d-single": "single"}
    for key, peers in sorted(by_cell.items(), key=str):
        if len(key) != 2 or key[0] not in names or not str(key[1])[-1].isdigit():
            continue
        w = int(str(key[1]).rsplit("-", 1)[-1])
        front = sorted((peers[i] for i in pareto(peers)), key=lambda e: (e["n"], -e["k"], -e["d"]))
        lines = [f"[[{e['n']},{e['k']},{e['d']}]] w={e['w']}" for e in front]
        (OUT / f"frontier_{names[key[0]]}_w{w}.txt").write_text(
            f"{key}: {len(peers)} codes, {len(front)} nondominated\n" + "\n".join(lines) + "\n", encoding="utf-8")
        print(key, len(peers), "codes,", len(front), "nondominated")


if __name__ == "__main__":
    main()
