# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""The [[264,8,12]] winner's tile on the paper's 10 x 10 bulk (the paper reports [[288,8,14]] as its best there)."""
import json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[4]))
from qec_search import gf2, tile
from qec_search.distqldpc import solve_many
if __name__ == "__main__":
    rows = [json.loads(l) for l in open(pathlib.Path(__file__).parents[1] / "proofs.jsonl")]
    r = next(r for r in rows if r.get("d_exact") == 12 and r["n"] == 264)
    hx, hz, _ = tile.build([tuple(e) for e in r["tile"]], 3, 10, 10)
    sol = solve_many({"paper10x10_264tile": (hx.astype("int8"), hz.astype("int8"))}, 3600, 8)["paper10x10_264tile"]
    print(r["key"], hx.shape[1], hx.shape[1] - gf2.rank(hx) - gf2.rank(hz), sol["sides"], flush=True)
