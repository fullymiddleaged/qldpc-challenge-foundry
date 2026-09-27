# CLAUDE.md — qldpc-lab

Project memory for Claude Code. Read `HANDOFF.md` for the full story; this file is the short, always-on version.

## What this repo is
Our own research code for **quantum LDPC codes** (bivariate bicycle / BB family first): search for codes, lay them out on a 2D chip, certify their distance, simulate them, and submit good ones to the Unitary Foundation **qLDPC Challenge** leaderboard. It lives in `lab/` on the main branch of our fork of the challenge repo (`origin` = fullymiddleaged/qldpc-challenge-foundry, `upstream` = unitaryfoundation/qldpc-challenge). The fork's `main` = upstream's files + `lab/`: never edit upstream's files on main, so `git merge upstream/main` stays conflict-free. Commands below run from `lab/`.

## Commands
- Setup (once): `bash scripts/setup_workspace.sh`
- Activate env: `source .venv/bin/activate` (Windows Git Bash: `source .venv/Scripts/activate`)
- Tests: `pytest -q` (no pytest? `python tests/run_all.py`) — about 1 minute, must pass before every commit
- Self-test incl. Stim circuits: `python -m qec_search.selftest`
- Search: `python -m qec_search.search --out runs/<name> --hours 2 [--mixed] [--terms 3,4] [--seed-known]`
- Report: `python -m qec_search.summarize runs/<name>`
- Export best codes + bilayer layouts: `python -m qec_search.export runs/<name> --layout --dest submissions/`
- Prove exact distance: `python -m qec_search.exact_distance --genome '<json>' --claimed <d> --hours <h> --out results/proofs/<file>.json`
- SAT certificate (upstream's encoding, resumable, from the repo root): `uv run --with pycryptosat --with python-sat python lab/scripts/certify_upstream.py <npz...> --d <d>`; progress: add `--status`
- Leaderboard frontier: `(cd .. && uv run python cli/qldpc.py targets --n 300 --top 10)` (or `./qldpc` from Git Bash; set `PYTHONUTF8=1` on Windows)
- Sync upstream: `git fetch upstream && git merge upstream/main` (the GitHub "Sync fork" button no longer works; Actions are disabled on the fork)
- Submit (dry run → real): `python scripts/submit.py <candidate> --handle @<user> [--note <md>] [--for-real]`. It runs upstream's own `qldpc submit` in a worktree `.worktrees/<candidate>` detached at `upstream/main`, so a PR never carries `lab/`. Never run `qldpc submit --open-pr` from the fork's main: it branches from the current HEAD.

## Layout
- `qec_search/` — the package (pure numpy/scipy core; `stim` + `ldpc` optional for simulation)
  - `gf2.py` GF(2) algebra · `bbcode.py` BB construction, distance upper bound, `components()`, `signature()` · `known_codes.py` IBM references
  - `decoders.py` ldpc BP-OSD / numpy fallback · `evaluate.py` Monte Carlo · `circuits.py` Stim memory circuits (**not yet validated on real Stim output beyond the self-test**)
  - `search.py` evolutionary search · `summarize.py` report · `verify.py` circuit-level slow loop
  - `layout.py` bilayer layout by simulated annealing · `export.py` submission files · `exact_distance.py` MILP proof
- `results/candidates/` submission-ready `.npz` (hx, hz, coords) + `manifest.json` with a **status** per code
- `results/proofs/` distance certificates · `results/logs/` checker + run logs · `colab/` Colab notebooks for heavy jobs (they clone the fork; see `colab/README.md`)
- `scripts/` setup + submit · `tests/` regression tests

## Rules (learned the hard way — do not skip)
1. **Distance claims:** anything from `distance_upper_bound` or the challenge's RIS search is an *upper bound*. Call a distance "exact" only with a finished MILP proof in `results/proofs/`. Our quick bound over-estimated twice at n=288 (claimed 24 → 22, 20 → 18).
2. **Duplicates:** the canonical key misses some equivalences. Use `code.signature()` (spectrum fingerprint) to merge relabelled copies; use `code.components() > 1` to drop codes that are several copies of a smaller code (e.g. [[288,24,12]] = 2 × gross).
3. **Attribution:** never submit other people's codes (IBM's gross code etc.) as records; a trivial layout of a known code is not a contribution. `submit.py` refuses `do_not_submit` / `hold` statuses.
4. **"Record" means record on the leaderboard under its rules**, not a physics result. Check the literature before calling anything new. Word construction notes modestly.
5. **Leaderboard rules:** interaction radius = max over checks of the largest pairwise distance between that check's *data* qubits; bilayer cap 7.0 (single layer 4.0); ≤ `layers` qubits per site; distinct sites ≥ 1 apart; n ≤ 700; ≤ 600 coordinates. Pareto dominance on (n↓, k↑, d↑, w↓) within a (locality, weight) cell. Non-records are accepted.
6. **One JSON per PR**, each on its own branch cut from `upstream/main` (`submit.py` does this). Follow upstream `AGENTS.md` for the PR body and research note: every path cited must exist in the PR's own tree, so `lab/...` paths can't be cited there.
7. **Long jobs** (MILP proofs, multi-hour searches): run in the background with a log file (`nohup … > log 2>&1 &`) or on Colab; never block the session on them.
8. Before claiming any number in docs or PRs, re-run the check that produces it.

## Style
Pure functions, numpy arrays of uint8 for GF(2) matrices, `from __future__ import annotations`, small CLIs with argparse. Keep results reproducible: every result file records genome, parameters, method, and date.
