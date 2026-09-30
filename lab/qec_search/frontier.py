# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""A board cell's Pareto frontier and the screening every hunt does against it. Every hunt imports this; add screening
speed-ups here, not in a script.

The frontier files are what scripts/dump_board_cells.py writes (results/hunt/frontier_<single|bilayer>_w<w>.txt,
lines '[[n,k,d]] w=W'). A code (n, k, d, w) is dominated when some entry is at least as good on every axis (n lower,
k and d higher, w lower); equal parameters count as dominated, because a tie is no record.

screen_distance is the one screen: the frontier fixes the smallest d a code with this (n, k, w) would need (bar), so
the randomised distance search stops at the first logical lighter than that. On 30 Sep a tile screen without it spent
about an hour on 70 tiles, nearly all rejected.
"""
from __future__ import annotations

import pathlib
import re

HUNT = pathlib.Path(__file__).resolve().parents[1] / "results" / "hunt"

Entry = tuple[int, int, int, int]


def parse(text: str) -> list[Entry]:
    """(n, k, d, w) for every '[[n,k,d]] w=W' entry in a frontier file or in qldpc targets output."""
    return [tuple(map(int, m)) for m in re.findall(r"\[\[(\d+),(\d+),(\d+)\]\]\s+w=(\d+)", text)]


def load(cell: str, w: int) -> list[Entry]:
    """cell is 'single' or 'bilayer'."""
    return parse((HUNT / f"frontier_{cell}_w{w}.txt").read_text(encoding="utf-8"))


def dominated(c: Entry, frontier: list[Entry]) -> bool:
    n, k, d, w = c
    return any(fn <= n and fk >= k and fd >= d and fw <= w for fn, fk, fd, fw in frontier)


def bar(n: int, k: int, w: int, frontier: list[Entry]) -> int:
    """Smallest d at which (n, k, d, w) is not dominated."""
    return 1 + max((fd for fn, fk, fd, fw in frontier if fn <= n and fk >= k and fw <= w), default=0)


def screen_distance(code, w: int, frontier: list[Entry], trials: int) -> int | None:
    """Upper bound on the code's distance, or None as soon as a logical below the frontier's bar is found.

    A cheap pass (trials / 10) runs first; the full pass runs only for codes that survive it. Both stop early at
    bar - 1, since a lighter logical already settles that the code is dominated."""
    n, k = code.n, code.k
    need = bar(n, k, w, frontier)
    if need > n:
        return None
    d = code.distance_upper_bound(trials=max(5, trials // 10), stop_at=need - 1)
    if d < need:
        return None
    d = min(d, code.distance_upper_bound(trials=trials, seed=1, stop_at=need - 1))
    return None if d < need else d
