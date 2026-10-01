# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Is a winning tile the paper's? arXiv:2504.09171's best 3x3 weight-8 tile gives [[288,8,14]] on a 10x10 bulk."""
import json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[4]))
from qec_search import gf2, tile
from qec_search.distqldpc import solve_many

if __name__ == "__main__":
    winners = {"e603e2da3b42": None, "a6374124aad4": None, "e8f30d7c1315": None, "67d6bb6dcbd5": None}
    rows = [json.loads(l) for l in open(pathlib.Path(__file__).parents[1] / "proofs.jsonl")]
    codes, keys = {}, {}
    for r in rows:
        stem = pathlib.Path(r["npz"]).stem
        h = next((w for w in winners if w in stem), None)
        if h and r.get("d_exact"):
            hx, hz, _ = tile.build([tuple(e) for e in r["tile"]], 3, 10, 10)
            codes[f"paper10x10_{h}"] = (hx.astype("int8"), hz.astype("int8"))
            keys[f"paper10x10_{h}"] = r["key"]
    sol = solve_many(codes, 3600, 8)
    out = {}
    for name, (hx, hz) in codes.items():
        k = hx.shape[1] - gf2.rank(hx) - gf2.rank(hz)
        out[name] = {"key": keys[name], "n": hx.shape[1], "k": k, "sides": sol[name]["sides"]}
        print(name, keys[name], hx.shape[1], k, {s: v for s, v in sol[name]["sides"].items()}, flush=True)
    (pathlib.Path(__file__).parent / "result.json").write_text(json.dumps(out, indent=1))
