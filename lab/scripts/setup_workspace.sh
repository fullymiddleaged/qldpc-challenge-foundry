#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
# One-time setup for the lab inside our fork of unitaryfoundation/qldpc-challenge.
#
#   usage:  bash lab/scripts/setup_workspace.sh        (from anywhere in the fork)
#
# Layout: the fork's main = upstream's files + our lab/ folder. Submissions never branch
# from main; lab/scripts/submit.py cuts them from upstream/main in .worktrees/<name>.
#
# Works on macOS, Linux, WSL and Git Bash on Windows.
set -euo pipefail

LAB="$(cd "$(dirname "$0")/.." && pwd)"
FORK="$(dirname "$LAB")"

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }

say "1/5  upstream remote and local excludes"
cd "$FORK"
git remote get-url upstream >/dev/null 2>&1 || \
  git remote add upstream https://github.com/unitaryfoundation/qldpc-challenge.git
git fetch upstream --quiet
grep -qx '.worktrees/' .git/info/exclude 2>/dev/null || echo '.worktrees/' >> .git/info/exclude
git remote -v | sed 's/^/  /'

say "2/5  Python environment for the lab (lab/.venv, includes uv for the challenge CLI)"
cd "$LAB"
PY="${PYTHON:-python3}"
if [ ! -d .venv ]; then "$PY" -m venv .venv; fi
if [ -f .venv/bin/activate ]; then . .venv/bin/activate; else . .venv/Scripts/activate; fi
python -m pip install --quiet --upgrade pip
python -m pip install --quiet -e ".[sim,dev]" uv

say "3/5  Self-test (IBM codes, decoder, Stim circuits)"
python -m qec_search.selftest

say "4/5  Unit tests"
python -m pytest -q

say "5/5  Challenge tool smoke test (current frontier, n <= 300)"
( cd "$FORK" && PYTHONUTF8=1 uv run --quiet python cli/qldpc.py targets --n 300 --top 3 ) || \
  echo "  (challenge tool did not run; see lab/HANDOFF.md troubleshooting)"

say "Done."
