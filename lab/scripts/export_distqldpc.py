# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Write one CSS side of a code as DistQLDPC input, for a manual run. The encoding (and every speed-up: duality,
orbital cubes) lives in qec_search/distqldpc.py; use qec_search.distqldpc.solve_many to solve rather than this script.

  Hx := h_opp        z must lie in ker(h_opp)
  Hz := I_n          forces every x variable to 0
  Gz := pairing_set  z must anticommute with some logical of the other type (all k classes, checked)
  Gx := one zero row (a zero row forces its selector off; the loader rejects an empty file)
--cube takes a certify_sym cube (1-based literals): negative literals become Hx rows e_i, positive ones are printed
as the -one-z list for the patched solver (scripts/distqldpc-one-z.patch).

    python scripts/export_distqldpc.py ../codes/180-18-14.json results/distqldpc --cube 91,-1,-2 --tag r1
"""
from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from qec_search.certify_sym import load_code  # noqa: E402
from qec_search.distqldpc import one_sided, write  # noqa: E402,F401  (one_sided/write re-exported for callers)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("code", type=pathlib.Path, help="board JSON or candidate .npz")
    ap.add_argument("out_dir", type=pathlib.Path)
    ap.add_argument("--side", choices=["X", "Z"], default="X")
    ap.add_argument("--cube", default="", help="comma-separated 1-based literals from a certify_sym status file")
    ap.add_argument("--tag", default="", help="suffix for the file stem")
    args = ap.parse_args()
    lits = [int(t) for t in args.cube.split(",") if t.strip()]
    hx, hz, _ = load_code(args.code)
    prefix = args.out_dir / f"{args.code.stem}_{args.side}{'_' + args.tag if args.tag else ''}"
    write(one_sided(hx, hz, args.side, [-v - 1 for v in lits if v < 0]), prefix)
    ones = [v - 1 for v in lits if v > 0]
    print(prefix, f"-one-z={','.join(map(str, ones))}" if ones else "")


if __name__ == "__main__":
    main()
