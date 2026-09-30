# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Audit the board's SAT distance certificates (certs/*.json made by upstream verify/sat_certify.py).

sat_certify pairs candidate logicals against sat_certify._logicals(H_opp, H_same), which takes rows
rank(H_same):rank(R) of the RREF of [H_same; ker H_opp]. Those rows need not complement rowspace(H_same). When
they span fewer than k logical classes, a nontrivial logical that commutes with every row escapes the query, so
UNSAT no longer proves the distance. Historical: upstream fixed _logicals in #2343 (28 Sep 2026, our issue #2273),
so against current upstream Stage 1 finds no deficient pairing sets; results/logs/sat_cert_audit.json records the
audit as run against the old code.

Stage 1 (static, seconds): for each SAT certificate, count the logical classes its pairing set spans on each side.
Stage 2 (--recertify): re-prove each deficient certificate with lab/qec_search/certify_sym.py. That uses a complete
logical basis, checked to span k classes. The outcomes: PROVEN (the board's distance holds), REFUTED (a lighter
logical exists, with its witness), or "in progress" (the budget ran out).

Run from the repo root:
    python lab/scripts/audit_sat_certs.py                                 # stage 1 only
    python lab/scripts/audit_sat_certs.py --recertify --max-cpu-hours 0.5  # stage 2, smallest codes first
Output: lab/results/logs/sat_cert_audit.json (rewritten as it goes; stage 2 resumes from certify_sym status files)
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lab"))
sys.path.insert(0, str(ROOT / "verify"))

import gf2  # noqa: E402  (upstream's verify/gf2.py)
import sat_certify  # noqa: E402  (upstream's verify/sat_certify.py, the code under audit)

from qec_search import certify_sym  # noqa: E402

OUT = ROOT / "lab" / "results" / "logs" / "sat_cert_audit.json"


def spanned_classes(doc: dict) -> dict:
    """k, and per side the classes spanned by the pairing set sat_certify.certify() actually used."""
    n = doc["n"]
    hx, hz = sat_certify._matrix(doc["checks"]["X"], n), sat_certify._matrix(doc["checks"]["Z"], n)
    k = n - gf2.rank(hx) - gf2.rank(hz)
    sides = {}
    for side, h_same, h_opp in (("X", hx, hz), ("Z", hz, hx)):
        t = sat_certify._logicals(h_opp, h_same)            # the exact call in sat_certify.certify()
        sides[side] = int(gf2.rank(np.vstack([h_opp, t])) - gf2.rank(h_opp)) if len(t) else 0
    return {"k": int(k), "spanned": sides}


def stage1() -> list[dict]:
    rows = []
    for cert in sorted((ROOT / "certs").glob("*.json")):
        c = json.loads(cert.read_text())
        code = ROOT / "codes" / cert.name
        if "SAT" not in str(c.get("solver", "")) or not code.exists():
            continue
        doc = json.loads(code.read_text())
        s = spanned_classes(doc)
        rows.append({"code": cert.name, "n": doc["n"], "k": s["k"], "d": doc["distance"]["d"],
                     "family": doc.get("family"), "cert_solver": c.get("solver"), "cert_d_exact": c.get("d_exact"),
                     "spanned": s["spanned"], "deficient": any(v < s["k"] for v in s["spanned"].values())})
    return rows


def save(rows: list[dict]) -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    bad = [r for r in rows if r["deficient"]]
    summary = {"sat_certs": len(rows), "deficient": len(bad),
               "recertified": {v: sum(r.get("recertify", {}).get("verdict") == v for r in bad)
                               for v in ("PROVEN", "REFUTED", "in progress")}}
    OUT.write_text(json.dumps({"updated": dt.datetime.now().isoformat(timespec="seconds"), "summary": summary,
                               "codes": rows}, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recertify", action="store_true", help="stage 2: re-prove deficient certificates")
    ap.add_argument("--max-cpu-hours", type=float, default=0.5, help="solver budget per code in stage 2")
    ap.add_argument("--tlim", type=float, default=120, help="seconds per cube before it is split")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--only", nargs="*", help="restrict stage 2 to these code files (e.g. 140-6-14.json)")
    args = ap.parse_args()

    rows = stage1()
    save(rows)
    bad = [r for r in rows if r["deficient"]]
    print(f"{len(rows)} SAT certificates, {len(bad)} with a pairing set spanning fewer than k classes", flush=True)
    if not args.recertify:
        return
    todo = sorted((r for r in bad if not args.only or r["code"] in args.only), key=lambda r: (r["n"], r["code"]))
    for r in todo:
        code = ROOT / "codes" / r["code"]
        st = certify_sym.certify([code], None, 2, True, args.tlim, args.workers, tag="_audit", quiet=True,
                                 max_cpu_hours=args.max_cpu_hours)[code]
        summ = certify_sym.summarize(st)
        r["recertify"] = {"verdict": summ["verdict"], "cpu_hours": summ["cpu_hours"],
                          "status_file": str(certify_sym.status_path(code, r["d"], "_audit").relative_to(ROOT))}
        wit = [c["witness"] for c in st["cubes"].values() if c.get("status") == "SAT"]
        if wit:
            r["recertify"]["lighter_logical"] = {"weight": len(wit[0]), "support": wit[0]}
        save(rows)
        print(f"{r['code']:>24} n={r['n']} k={r['k']} d={r['d']} spanned={r['spanned']} -> {summ['verdict']} "
              f"({summ['cpu_hours']} cpu-h)" + (f" LIGHTER LOGICAL weight {len(wit[0])}" if wit else ""), flush=True)


if __name__ == "__main__":
    main()
