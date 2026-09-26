# qldpc-lab

Search, lay out, certify and submit **quantum LDPC codes**, starting with the bivariate bicycle family behind IBM's gross code. The work targets both the research frontier and the [Unitary Foundation qLDPC Challenge](https://github.com/unitaryfoundation/qldpc-challenge) leaderboard.

```
search ──► dedupe ──► 2D layout ──► exact distance proof ──► challenge checker ──► PR
(search.py)  (signature,  (layout.py)   (exact_distance.py)     (upstream verifier)    (scripts/submit.py)
             components)
```

This folder lives on the `main` branch of our fork of the challenge; everything outside `lab/` is upstream's.

**Start here:** `HANDOFF.md` has the full context, findings and roadmap. `CLAUDE.md` holds the working rules for Claude Code.

## Quick start
```bash
bash scripts/setup_workspace.sh   # upstream remote, lab/.venv, self-test, tests
source .venv/bin/activate        # Git Bash on Windows: source .venv/Scripts/activate
python -m qec_search.search --out runs/mixed1 --hours 3 --mixed --seed-known
python -m qec_search.summarize runs/mixed1
python -m qec_search.export runs/mixed1 --layout --dest submissions/
python scripts/submit.py <candidate> --handle @you             # dry run through the official checker, on upstream/main
```

## Current results (2026-09-26)
| Code | Cell | Distance | Status |
|---|---|---|---|
| [[216,8,16]] | 2D-local bilayer, w ≤ 6 | proven exact | ready |
| [[288,12,16]] | 2D-local bilayer, w ≤ 6 | proven exact | ready |
| [[288,8,20]] | 2D-local bilayer, w ≤ 6 | ≤ 20 (proof pending) | hold |

Heavy runs: `colab/`. Files: `results/candidates/`, proofs: `results/proofs/`, checker logs: `results/logs/`.
