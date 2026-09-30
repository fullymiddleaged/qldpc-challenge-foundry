# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""qec_search.frontier: the dominance test, the bar a new code must clear, and the early-stopping screen."""
from __future__ import annotations

from qec_search import frontier as fr
from qec_search.bbcode import BBGenome, build_bb

GROSS = build_bb(BBGenome(12, 6, ((3, 0), (0, 1), (0, 2)), ((0, 3), (1, 0), (2, 0))))    # [[144,12,12]]
F = fr.parse("[[112,12,12]] w=8\n[[192,12,12]] w=6\n[[144,8,12]] w=6\n")


def test_bar_is_one_above_the_best_comparable_entry():
    assert fr.bar(144, 12, 6, F) == 1          # nothing with n <= 144, k >= 12 and w <= 6
    assert fr.bar(144, 12, 8, F) == 13         # [[112,12,12]] w=8 must be beaten
    assert fr.bar(200, 8, 6, F) == 13          # both w=6 entries cover it at 12
    assert fr.bar(144, 8, 6, F) == 13          # a tie is dominated, so 12 is not enough


def test_bar_agrees_with_dominated_at_the_boundary():
    for n, k, w in [(144, 12, 8), (200, 8, 6), (144, 8, 6)]:
        b = fr.bar(n, k, w, F)
        assert fr.dominated((n, k, b - 1, w), F) and not fr.dominated((n, k, b, w), F)


def test_screen_keeps_the_gross_code_where_it_is_new():
    assert fr.screen_distance(GROSS, 6, F, trials=30) == 12


def test_screen_rejects_and_stops_early_where_it_is_beaten(monkeypatch):
    calls = []
    real = type(GROSS).distance_upper_bound

    def spy(self, trials=300, seed=0, which="both", stop_at=None):
        calls.append(stop_at)
        return real(self, trials=trials, seed=seed, which=which, stop_at=stop_at)
    monkeypatch.setattr(type(GROSS), "distance_upper_bound", spy)
    assert fr.screen_distance(GROSS, 8, F, trials=30) is None       # [[112,12,12]] w=8 dominates it
    assert calls == [12]                        # one cheap pass, told to stop at bar - 1 = 12, and no full pass
