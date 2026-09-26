# CLAUDE.md: our fork of the qLDPC Challenge

- `lab/` is our research project; read `lab/CLAUDE.md` before working in it. Everything else here is upstream's
  (unitaryfoundation/qldpc-challenge); follow `AGENTS.md` for how the board works.
- Don't edit upstream files on `main`; keep our changes under `lab/` so `git merge upstream/main` stays clean.
- Leaderboard submissions go through `lab/scripts/submit.py`, which branches from `upstream/main` in
  `.worktrees/`. Never run `qldpc submit --open-pr` from `main`, because it would ship `lab/` in the PR.
- Submit or push only when the user says so.
