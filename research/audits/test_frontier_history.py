"""Tests for the frontier replay (``frontier_history.py``).

The script exists so that other people can quote its numbers, and its numbers
are only as good as its replay of ``codes/``, so the invariant that matters is
stated once and asserted directly: whatever walk ``replay()`` takes through
history, it must end where the repository is. The final state is compared to
``git ls-tree -r HEAD -- codes/`` read as ``(n, k, d)`` triples, CSS only --
three seconds of git, and the exact check that would have caught a replay that
applies commits in an order the ancestry never had.

The rest pins what is wrong silently if it regresses: a rename is read at its
destination (git reports a rewrite of an unrelated file as a rename when the
similarity check matches, so carrying the source's value over invents one no
blob has); the three-axis frontier is a strict antichain; the hypervolume is
the union volume, not the sum; a landing is dated by when it landed; and the
region between two staircases is split into gains and losses correctly.
"""

import math
import os
import shutil
import subprocess
import sys
import types

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import frontier_history as fh  # noqa: E402


def head_css_state(repo):
    """``{path: (n, k, d)}`` for every CSS code in ``git ls-tree -r HEAD``."""
    out = subprocess.run(["git", "ls-tree", "-r", "HEAD", "--", "codes/"],
                         cwd=repo, capture_output=True, text=True,
                         check=True).stdout
    paths = [line.split("\t")[1] for line in out.splitlines()
             if line.split("\t")[1].endswith(".json")]
    blobs = fh.fetch_params([("HEAD", p) for p in paths], repo)
    return {path: value[:3] for (_rev, path), value in blobs.items()
            if value and value[3] == "CSS"}


def test_replay_ends_where_head_is():
    """The state the walk reaches is the state the repository is in.

    This is the whole trust argument for the charts: a replay that ends at a
    board HEAD never had -- one file too many, one file too few, or the same
    file carrying another commit's ``(n, k, d)`` -- reports a frontier that
    describes the walk rather than the board.
    """
    rows, state = fh.replay(fh._REPO)
    assert rows, "no board-moving commit found"
    assert state == head_css_state(fh._REPO)
    # The last row is the board at the last commit that moved it, and nothing
    # after it moved the board, so the row and HEAD agree on size too.
    assert rows[-1][2] == len(state)
    assert set(rows[-1][5]) == set(fh.frontier(set(state.values())))
    # Every blob the walk asked for parsed: a silent skip would drop a code
    # from the replay without dropping it from HEAD.
    assert fh.replay.unparsed == 0


def test_rename_is_read_at_its_destination(monkeypatch):
    """A rename takes the destination's own blob, not the source's value."""
    log = "\n".join([
        "C aaaa1111 2026-01-01T00:00:00+00:00",
        "M\tcodes/old.json",
        "C bbbb2222 2026-01-02T00:00:00+00:00",
        "R100\tcodes/old.json\tcodes/new.json",
    ])
    seen = {}

    def fake_git(args, cwd):
        return types.SimpleNamespace(returncode=0, stdout=log.encode(),
                                     stderr=b"")

    def fake_fetch(requests, cwd):
        seen["requests"] = list(requests)
        return {
            ("aaaa1111", "codes/old.json"): (64, 8, 4, "CSS"),
            ("bbbb2222", "codes/new.json"): (100, 8, 6, "CSS"),
        }

    monkeypatch.setattr(fh, "git", fake_git)
    monkeypatch.setattr(fh, "fetch_params", fake_fetch)

    rows, state = fh.replay("unused")
    assert ("bbbb2222", "codes/new.json") in seen["requests"]
    assert state == {"codes/new.json": (100, 8, 6)}
    assert [r[2] for r in rows] == [1, 1]
    # The rename moved (64, 8, 4) off the frontier and the row names it.
    assert rows[1][7] == "-[[64,8,4]]"


def test_a_distance_revision_is_named_as_one(monkeypatch):
    """An in-place d correction reads ``[[n,k,d]]->d'`` in the row's note."""
    log = "\n".join([
        "C aaaa1111 2026-01-01T00:00:00+00:00",
        "A\tcodes/a.json",
        "C bbbb2222 2026-01-02T00:00:00+00:00",
        "M\tcodes/a.json",
    ])

    def fake_git(args, cwd):
        return types.SimpleNamespace(returncode=0, stdout=log.encode(),
                                     stderr=b"")

    def fake_fetch(requests, cwd):
        return {
            ("aaaa1111", "codes/a.json"): (882, 18, 30, "CSS"),
            ("bbbb2222", "codes/a.json"): (882, 18, 29, "CSS"),
        }

    monkeypatch.setattr(fh, "git", fake_git)
    monkeypatch.setattr(fh, "fetch_params", fake_fetch)
    rows, _ = fh.replay("unused")
    assert rows[1][7] == "[[882,18,30]]->29"
    assert rows[1][6] < rows[0][6]


def test_params_reads_d_from_any_schema_shape():
    assert fh.params(b'{"n": 7, "k": 1, "distance": {"d": 3}}') == \
        (7, 1, 3, "CSS")
    assert fh.params(b'{"n": 7, "k": 1, "distance": {"X": {"value": 3}, '
                     b'"Z": {"value": 5}}}') == (7, 1, 3, "CSS")
    assert fh.params(b'{"n": 7, "k": 1, "distance": 3, '
                     b'"code_type": "stabilizer"}') == (7, 1, 3, "stabilizer")
    assert fh.params(b'{"n": 7, "k": 1}') is None
    assert fh.params(b"not json") is None


def test_frontier_is_a_strict_antichain():
    """No triple on the frontier is beaten on all three axes; every other is."""
    points = [(4, 4, 2), (4, 8, 2), (8, 6, 2), (8, 2, 4), (16, 10, 2),
              (32, 9, 2), (16, 10, 3), (4, 8, 1)]
    front = fh.frontier(points)
    assert front == [(4, 8, 2), (8, 2, 4), (16, 10, 3)]

    def beats(a, b):
        return (a[0] <= b[0] and a[1] >= b[1] and a[2] >= b[2]
                and a != b)

    for p in front:
        assert not any(beats(q, p) for q in points)
    for p in set(points) - set(front):
        assert any(beats(q, p) for q in points)


def test_frontier_matches_brute_force():
    import random
    rng = random.Random(7)
    points = {(rng.randrange(4, 300, 2), rng.randrange(1, 60),
               rng.randrange(2, 40)) for _ in range(400)}
    fast = set(fh.frontier(points))
    slow = {p for p in points
            if not any(q[0] <= p[0] and q[1] >= p[1] and q[2] >= p[2]
                       and q != p for q in points)}
    assert fast == slow


def test_hypervolume_is_the_union_not_the_sum():
    # One box: log2(1000/125) x (1 + log2 4) x (1 + log2 4) = 3 x 3 x 3.
    assert math.isclose(fh.hypervolume([(125, 4, 4)], 1000), 27.0)
    # A dominated point adds nothing; a disjoint corner adds its own volume
    # minus the overlap.
    assert math.isclose(fh.hypervolume([(125, 4, 4), (250, 2, 2)], 1000), 27.0)
    both = fh.hypervolume([(125, 4, 4), (125, 16, 1)], 1000)
    # Boxes 3x3x3 and 3x5x1 overlap in 3x3x1: 27 + 15 - 9.
    assert math.isclose(both, 33.0)
    # At or past the reference n a code contributes nothing.
    assert fh.hypervolume([(1000, 100, 100)], 1000) == 0.0
    assert fh.hypervolume([], 1000) == 0.0


def test_landing_is_dated_by_its_utc_day():
    assert fh.utc_day("2026-09-29T07:54:29-07:00") == "2026-09-29"
    assert fh.utc_day("2026-09-29T23:30:00-07:00") == "2026-09-30"
    assert fh.utc_day("2026-09-30T00:30:00+02:00") == "2026-09-29"


def test_floor_staircase_keeps_only_codes_clearing_the_floor():
    front = [(12, 2, 4), (16, 6, 4), (20, 8, 4), (24, 4, 8), (64, 18, 8),
             (64, 12, 12)]
    assert fh.floor_staircase(front, 4) == [(12, 2, 4), (16, 6, 4),
                                            (20, 8, 4), (64, 18, 8)]
    assert fh.floor_staircase(front, 8) == [(24, 4, 8), (64, 18, 8)]
    assert fh.floor_staircase(front, 12) == [(64, 12, 12)]
    assert fh.floor_staircase(front, 40) == []


def test_region_rects_split_gain_from_loss():
    old = [(10, 2, 4), (20, 6, 4), (40, 10, 4)]
    new = [(10, 4, 4), (20, 5, 4), (40, 10, 4), (80, 20, 4)]
    rects = fh.region_rects(old, new, 100)
    assert rects == [("gain", 10, 20, 2, 4), ("loss", 20, 40, 5, 6),
                     ("gain", 80, 100, 10, 20)]
    assert fh.region_rects(old, old, 100) == []


def test_snapshots_keep_the_last_row_of_each_bucket():
    rows = [("2026-01-01", "a", 10, 1, 1, ((4, 4, 2),), 1.0, ""),
            ("2026-01-02", "b", 11, 2, 2, ((4, 4, 2), (8, 6, 2)), 2.0, ""),
            ("2026-01-09", "c", 12, 2, 2, ((4, 4, 2), (8, 6, 2)), 2.0, "")]
    assert [r[0] for r in fh.snapshots(rows, "week")] == \
        ["2026-01-02", "2026-01-09"]
    assert [r[0] for r in fh.snapshots(rows, "day")] == [r[0] for r in rows]
    assert [r[0] for r in fh.snapshots(rows, "month")] == ["2026-01-09"]


@pytest.mark.parametrize("bucket", ["day", "week", "month"])
def test_bucket_key_is_monotone_within_a_bucket(bucket):
    """Dates inside one bucket share a key and sort with it, not against it."""
    keys = [fh.bucket_key(d, bucket) for d in
            ("2026-01-01", "2026-01-31", "2026-02-01")]
    assert keys == sorted(keys)


_ROWS = [("2026-01-05", "a", 2, 2, 2, ((12, 2, 4), (64, 12, 12)), 30.0, ""),
         ("2026-02-05", "b", 3, 3, 3, ((12, 2, 4), (24, 4, 8),
                                       (64, 12, 12)), 40.0, ""),
         ("2026-03-05", "c", 3, 2, 2, ((12, 2, 4), (24, 4, 8)), 35.0,
          "-[[64,12,12]]")]


def test_charts_render_from_synthetic_rows(tmp_path):
    snaps = fh.svg_frontier(_ROWS, tmp_path / "f.svg", "month", (4, 8))
    assert len(snaps) == 3
    svg = (tmp_path / "f.svg").read_text()
    assert "url(#loss)" in svg and "64,12,12" not in svg.split("<circle")[0]
    drops = fh.svg_history(_ROWS, tmp_path / "h.svg")
    assert len(drops) == 1
    assert "-[[64,12,12]]" in (tmp_path / "h.svg").read_text()


def test_frames_are_one_per_change_over_the_previous(tmp_path):
    # A fourth snapshot with the same staircases as the third folds into it.
    rows = _ROWS + [("2026-03-20", "d", 4, 2, 2, ((12, 2, 4), (24, 4, 8)),
                     35.0, "")]
    paths = fh.write_frames(rows, tmp_path / "frames", "day", (4, 8))
    assert [os.path.basename(p) for p in paths] == [
        "frame_0000_2026-01-05.svg", "frame_0001_2026-02-05.svg",
        "frame_0002_2026-03-05.svg"]
    first, second, third = (open(p, encoding="utf-8").read() for p in paths)
    assert "previous frame" not in first and "previous frame" in second
    assert "as of 2026-03-05" in third and "url(#loss)" in third
    assert "unchanged through 2026-03-20" in third and "4 CSS codes" in third
    # Axes are shared: the same tick labels appear in every frame.
    ticks = lambda s: [ln for ln in s.splitlines() if 'fill="#888"' in ln]  # noqa: E731
    assert ticks(first) == ticks(third)


def test_gif_has_one_frame_per_snapshot(tmp_path):
    Image = pytest.importorskip("PIL.Image")
    if shutil.which("rsvg-convert") is None:
        pytest.skip("rsvg-convert not on PATH")
    count = fh.write_gif(_ROWS, tmp_path / "f.gif", "day", (4, 8),
                         frame_ms=100, width=400)
    assert count == 3
    with Image.open(tmp_path / "f.gif") as im:
        assert im.format == "GIF" and im.n_frames == 3
