# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Turn a run directory into a compact report to send back to Claude.

    python -m qec_search.summarize runs/run1 [runs/run2 ...]

Writes <first run dir>/report_for_claude.md (human- and model-readable,
kept short enough to paste) and top_candidates.json (exact genomes, so any
result can be rebuilt and re-checked).
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from .known_codes import KNOWN, best_known_kd2n


def load(run_dirs):
    recs, metas = {}, []
    for d in run_dirs:
        d = Path(d)
        if (d / "meta.json").exists():
            metas.append(json.load((d / "meta.json").open()))
        if (d / "results.jsonl").exists():
            for line in (d / "results.jsonl").open():
                try:
                    r = json.loads(line)
                except Exception:
                    continue
                old = recs.get(r["key"])
                if old is None or (r.get("cc") and not old.get("cc")):
                    recs[r["key"]] = r
        if (d / "circuit.jsonl").exists():
            for line in (d / "circuit.jsonl").open():
                try:
                    c = json.loads(line)
                except Exception:
                    continue
                if c["key"] in recs:
                    recs[c["key"]].setdefault("circuit", []).append(c["result"])
    return list(recs.values()), metas


def fmt_ler(cc):
    if not cc:
        return "-"
    return f"{cc['ler']:.2e} ({cc['fails']}/{cc['shots']})"


def params(r):
    return f"[[{r['n']},{r['k']},{r.get('d_ub', '?')}]]"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args(argv)
    all_recs, metas = load(args.runs)
    n_split = sum(1 for r in all_recs if r.get("split"))
    n_dup = sum(1 for r in all_recs if r.get("dup_of"))
    variants = {}
    for r in all_recs:
        if r.get("dup_of"):
            variants[r["dup_of"]] = variants.get(r["dup_of"], 0) + 1
    recs = [r for r in all_recs if not r.get("split") and not r.get("dup_of")]
    for r in recs:
        r["variants"] = variants.get(r["key"], 0)
    out_dir = Path(args.runs[0])
    L = []
    w = L.append

    w("# QEC search report")
    w("")
    w("## Run setup")
    for m in metas:
        a = m.get("args", {})
        w(f"- started {m.get('started')} finished {m.get('finished', 'still running / interrupted')}; "
          f"decoder={m.get('decoder')} stim={m.get('stim')} ldpc={m.get('ldpc')}; "
          f"evals={m.get('evals_total', m.get('evals_this_run'))} (k=0: {m.get('zero_k_total', m.get('zero_k_this_run'))}); "
          f"terms={a.get('terms')} mixed={a.get('mixed')} ps={a.get('ps')} shots={a.get('shots')} "
          f"dist_trials={a.get('dist_trials')} seed={a.get('seed')}")
        tbs = m.get("tried_by_size")
        if tbs:
            w(f"- genomes tried per torus size: " + ", ".join(f"{k}:{v}" for k, v in sorted(tbs.items())))
    w(f"- codes with k>0 recorded: {len(all_recs)}; split into disjoint copies (dropped): {n_split}; "
      f"relabelled duplicates (merged): {n_dup}; distinct codes kept: {len(recs)}; with Monte Carlo: {sum(1 for r in recs if r.get('cc'))}")
    w("")

    # ----- sanity: known codes
    w("## Reproduction check (IBM reference codes found in this run)")
    kn = [r for r in recs if r.get("known")]
    if kn:
        for r in sorted(kn, key=lambda r: r["n"]):
            exp = next(p for name, _, p in KNOWN if name == r["known"])
            ok = "OK" if (r["n"], r["k"], r.get("d_ub")) == exp else f"MISMATCH expected {exp}"
            cc = r.get("cc") or [None]
            w(f"- {r['known']}: got {params(r)} {ok}; LER@p0 {fmt_ler(cc[0])}")
    else:
        w("- none of the reference codes were evaluated (use --seed-known to include them)")
    w("")

    # ----- frontier by n
    w("## Frontier: best k*d^2/n at each n (distance is an upper bound)")
    w("| n | best [[n,k,d]] | k*d^2/n | best known <= n | ratio | LER@p0 | genome |")
    w("|---|---|---|---|---|---|---|")
    by_n = defaultdict(list)
    for r in recs:
        if "kd2n" in r:
            by_n[r["n"]].append(r)
    for n in sorted(by_n):
        best = max(by_n[n], key=lambda r: (r["kd2n"], -(r.get("cc") or [{"ler": 1}])[0]["ler"]))
        ref = best_known_kd2n(n)
        ratio = f"{best['kd2n'] / ref:.2f}" if ref else "-"
        cc = (best.get("cc") or [None])[0]
        tag = f" ({best['known']})" if best.get("known") else ""
        w(f"| {n} | {params(best)}{tag} | {best['kd2n']:.2f} | {ref:.2f} | {ratio} | {fmt_ler(cc)} | `{best['pretty']}` |")
    w("")

    # ----- top list
    w(f"## Top {args.top} by k*d^2/n (with Monte Carlo)")
    w("| # | [[n,k,d]] | k*d^2/n | LER@p0 | LER@p1 | genome | note |")
    w("|---|---|---|---|---|---|---|")
    scored = [r for r in recs if r.get("cc")]
    scored.sort(key=lambda r: (r["kd2n"], -r["cc"][0]["ler"]), reverse=True)
    for i, r in enumerate(scored[: args.top], 1):
        cc = r["cc"]
        p1 = fmt_ler(cc[1]) if len(cc) > 1 else "-"
        note = r.get("known") or ("beats reference at this n" if r["kd2n"] > best_known_kd2n(r["n"]) + 1e-9 else "")
        if r.get("variants"):
            note = (note + f" (+{r['variants']} relabelled copies)").strip()
        w(f"| {i} | {params(r)} | {r['kd2n']:.2f} | {fmt_ler(cc[0])} | {p1} | `{r['pretty']}` | {note} |")
    w("")

    # ----- decoding anomalies: same (n,k,d), very different LER
    w("## Decoding anomalies (same [[n,k,d]], LER differs by >3x)")
    groups = defaultdict(list)
    for r in scored:
        groups[(r["n"], r["k"], r["d_ub"])].append(r)
    anyan = False
    for key, g in sorted(groups.items(), key=lambda kv: -kv[0][0]):
        g = [r for r in g if r["cc"][0]["fails"] >= 10]
        if len(g) < 2:
            continue
        lo = min(g, key=lambda r: r["cc"][0]["ler"])
        hi = max(g, key=lambda r: r["cc"][0]["ler"])
        if hi["cc"][0]["ler"] > 3 * lo["cc"][0]["ler"]:
            anyan = True
            w(f"- [[{key[0]},{key[1]},{key[2]}]]: best {fmt_ler(lo['cc'][0])} `{lo['pretty']}` vs worst {fmt_ler(hi['cc'][0])} `{hi['pretty']}`")
    if not anyan:
        w("- none")
    w("")

    # ----- structure stats
    w("## Structure statistics")
    k_by_size = defaultdict(list)
    for r in recs:
        k_by_size[f"{r['genome']['l']}x{r['genome']['m']}"].append(r)
    w("| torus | codes k>0 | k values seen | max d | max k*d^2/n |")
    w("|---|---|---|---|---|")
    for s in sorted(k_by_size, key=lambda s: (int(s.split('x')[0]) * int(s.split('x')[1]), s)):
        g = k_by_size[s]
        ks = sorted({r["k"] for r in g})
        w(f"| {s} | {len(g)} | {ks[:10]}{'...' if len(ks) > 10 else ''} | {max(r.get('d_ub', 0) for r in g)} | "
          f"{max(r.get('kd2n', 0) for r in g):.2f} |")
    w("")

    circ = [r for r in recs if r.get("circuit")]
    if circ:
        w("## Circuit-level results")
        w("| [[n,k,d]] | p | rounds | LER (block) | per round | detectors | genome |")
        w("|---|---|---|---|---|---|---|")
        for r in sorted(circ, key=lambda r: r["n"]):
            for c in r["circuit"]:
                w(f"| {params(r)} | {c['p']} | {c['rounds']} | {c['ler']:.2e} ({c['fails']}/{c['shots']}) | "
                  f"{c['ler_per_round']:.2e} | {c['num_detectors']} | `{r['pretty']}` |")
        w("")

    w("## Caveats")
    w("- d is an information-set upper bound, not a proof; verify finalists with an exact method.")
    w("- code-capacity numbers ignore measurement errors; only circuit-level results are decisive.")
    w("- 'beats reference' compares against 5 IBM codes only. The live frontier is the Unitary Foundation qLDPC Challenge "
      "leaderboard (unitaryfoundation.github.io/qldpc-challenge), which is far ahead of those 5 at weight 6.")

    text = "\n".join(L) + "\n"
    (out_dir / "report_for_claude.md").write_text(text)
    top = [{"key": r["key"], "params": params(r), "kd2n": r.get("kd2n"), "genome": r["genome"],
            "cc": r.get("cc"), "circuit": r.get("circuit"), "known": r.get("known")} for r in scored[:100]]
    (out_dir / "top_candidates.json").write_text(json.dumps(top, indent=1))
    print(text)
    print(f"wrote {out_dir / 'report_for_claude.md'} and {out_dir / 'top_candidates.json'}")


if __name__ == "__main__":
    main()
