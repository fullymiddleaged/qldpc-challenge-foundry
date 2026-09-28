# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Write a code as DistQLDPC input (github.com/guluchen/DistQLDPC, the MaxCDCL solver of arXiv:2606.12445).

DistQLDPC minimises Pauli weight over S^perp \\ S with both x and z variables. For a CSS code we pose one side only:
  Hx := h_opp        z must lie in ker(h_opp)
  Hz := I_n          forces every x variable to 0
  Gz := pairing_set  z must anticommute with some logical of the other type (all k classes, checked)
  Gx := one zero row (a zero row forces its selector off; the loader rejects an empty file)
so the optimum is exactly the side's distance, as in certify_sym. Side X means an X-type logical: x in ker(Hz).
--cube takes a certify_sym cube (1-based literals, e.g. "91,-1,-2"): each negative literal becomes an Hx row e_i, and
the positive ones are printed as the -one-z list for the patched solver (scripts/distqldpc-one-z.patch). The optimum is
then the lightest logical in that cube, so a cube is closed when it reaches the claimed d.

    python scripts/export_distqldpc.py ../codes/180-18-14.json results/distqldpc      # writes 180-18-14_X_{Hx,..}.txt
    python scripts/export_distqldpc.py ../codes/180-18-14.json results/distqldpc --cube 91,-1,-2 --tag r1
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from qec_search.certify_sym import load_code, pairing_set  # noqa: E402


def one_sided(hx: np.ndarray, hz: np.ndarray, side: str, zeros: list[int] = ()) -> dict[str, np.ndarray]:
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
