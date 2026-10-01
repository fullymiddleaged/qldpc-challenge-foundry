# CLAUDE.md — qldpc-lab

Project memory for Claude Code. Read `HANDOFF.md` for the full story; this file is the short, always-on version.

## What this repo is
**Aim:** Searching for quantum error-correcting codes that are both strong and buildable: laid out for real chips, with proven distances, and submitted to the qLDPC Challenge.

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
- Exact distance for any CSS code (orbital branching on Tanner-graph automorphisms, complete logical basis, resumable; needs `pip install -e .[cert]`): `python -m qec_search.certify_sym ../codes/<name>.json [--d <d>] [--depth 2] [--max-cpu-hours h]`; progress: `--status`. Supersedes `scripts/certify_upstream.py` (BB tori only)
- Fast exact distance: `qec_search.distqldpc.solve_many({name: (hx, hz)}, secs)` (DistQLDPC MaxSAT, arXiv:2606.12445, built in WSL at `~/src/DistQLDPC` + `scripts/distqldpc-one-z.patch`). It is the ONE solver path: one CSS side, X/Z duality, orbital cubes. Add speed-ups there, never in a script (the audit's own copy lacked the orbital split and left every n >= 168, d >= 16 code open). Minutes at d <= 18; no proof file, so it is a solver claim like upstream's certs. Logs in `results/distqldpc/logs/`
- Hunt new 2D-local bilayer entries (screen → layout → DistQLDPC proof, each stage resumable): `python scripts/hunt_bilayer.py {screen,layout,prove} --out results/hunt/<run> ...`; overnight driver `bash scripts/overnight_hunt.sh`. Needs the cell frontiers from `scripts/dump_board_cells.py` (run it in `.worktrees/frontier` at upstream/main; loads the board once, ~6 min, instead of `qldpc targets` per query)
- Screening against a board cell goes through `qec_search.frontier` (parse, dominated, `bar`, `screen_distance`): the distance search stops at the first logical below the bar a new code must clear. Every hunt imports it; add screening speed-ups there, never in a script (a tile screen without it took ~1 h for 70 tiles). Single-layer tile hunt: `python scripts/hunt_tile.py {screen,prove} --out results/hunt/<run> --B 3 --w 6`; overnight `bash scripts/overnight_tile.sh`
- Audit uncertified board distances (no file in upstream `certs/`) with DistQLDPC; any overclaim is a revision submission (TRACKS.md): `python scripts/audit_board_distances.py --local-only --n-max 400 --secs 900`; one line per code in `results/logs/board_distance_audit.jsonl`, reruns skip audited codes. Hunt rounds skip every genome listed in `results/hunt/*/screened.txt`; what each round covered is in `notes/hunt-*.md`
- Audit the board's SAT certificates (upstream `sat_certify._logicals` spanned < k classes on 255/280; fixed upstream by #2343 after our issue #2273, so this is historical): `python scripts/audit_sat_certs.py [--recertify --max-cpu-hours 0.25]` from `lab/` or the root; results in `results/logs/sat_cert_audit.json`
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
  - `automorphisms.py` Tanner-graph automorphisms (individualisation-refinement, every perm verified) · `certify_sym.py` symmetry-reduced SAT certifier
- `results/candidates/` submission-ready `.npz` (hx, hz, coords) + `manifest.json` with a **status** per code
- `results/proofs/` distance certificates · `results/logs/` checker + run logs · `colab/` Colab notebooks for heavy jobs (they clone the fork; see `colab/README.md`)
- `scripts/` setup + submit · `tests/` regression tests

## Rules (learned the hard way — do not skip)
1. **Distance claims:** anything from `distance_upper_bound` or the challenge's RIS search is an *upper bound*. Call a distance "exact" only with a finished proof in `results/proofs/`: MILP (`qec_search.exact_distance`) for BB codes; for codes without the BB torus symmetry (tile codes), a finished DistQLDPC optimum (`qec_search.distqldpc`) on every side, or on one side plus a verified X/Z duality, with the method named in the proof (user decision, 1 Oct 2026). Our quick bound over-estimated twice at n=288 (claimed 24 → 22, 20 → 18).
2. **Duplicates:** the canonical key misses some equivalences. Use `code.signature()` (spectrum fingerprint) to merge relabelled copies; use `code.components() > 1` to drop codes that are several copies of a smaller code (e.g. [[288,24,12]] = 2 × gross).
3. **Attribution:** never submit other people's codes (IBM's gross code etc.) as records; a trivial layout of a known code is not a contribution. A non-trivial layout of a known code may go in only credited, `not_our_code`, `novelty: known_parameters`, with the user's OK (the gross-code bilayer layout, PR #2394). Don't extend another board contributor's active line either: e-eight's B = 4 weight-8 tile (codes/578-18-20 and relatives) at new sizes is theirs to finish; we use it only as a benchmark and search our own tiles. `submit.py` refuses `do_not_submit` / `hold` statuses.
4. **"Record" means record on the leaderboard under its rules**, not a physics result. Check the literature before calling anything new. Word construction notes modestly.
5. **Leaderboard rules:** interaction radius = max over checks of the largest pairwise distance between that check's *data* qubits; bilayer cap 7.0 (single layer 4.0); ≤ `layers` qubits per site; distinct sites ≥ 1 apart; n ≤ 700; ≤ 600 coordinates. Pareto dominance on (n↓, k↑, d↑, w↓) within a (locality, weight) cell. Non-records are accepted.
6. **One JSON per PR**, each on its own branch cut from `upstream/main` (`submit.py` does this). Follow upstream `AGENTS.md` for the PR body and research note: every path cited must exist in the PR's own tree, so `lab/...` paths can't be cited there.
7. **Long jobs** (MILP proofs, multi-hour searches): run in the background with a log file (`nohup … > log 2>&1 &`) or on Colab; never block the session on them.
8. Before claiming any number in docs or PRs, re-run the check that produces it.

## Style
Pure functions, numpy arrays of uint8 for GF(2) matrices, `from __future__ import annotations`, small CLIs with argparse. Keep results reproducible: every result file records genome, parameters, method, and date.
