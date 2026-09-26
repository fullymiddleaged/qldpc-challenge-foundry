> Detailed usage notes for the `qec_search` package, carried over from the original qec-search README. Run commands from `lab/`; setup is `bash scripts/setup_workspace.sh` (see `lab/CLAUDE.md`).

# qec-search

An evolutionary search for quantum error-correcting codes in the **bivariate bicycle (BB)** family. This is the family behind IBM's [[144,12,12]] "gross code" and its fault-tolerance roadmap. Everything runs on an ordinary computer.

The idea: you run thousands of cheap simulated experiments, then send the compact report back to Claude. Claude looks for patterns and proposes where to search next or what to change, and you run again.

```
 you: run search ──► results.jsonl ──► summarize ──► report_for_claude.md ──► Claude
  ▲                                                                            │
  └──────────── new search settings / code changes / ideas ◄───────────────────┘
```

## What it does

For each candidate code (a torus size l×m plus two polynomials A, B):

1. **Builds** the check matrices H_X = [A|B], H_Z = [Bᵀ|Aᵀ] and computes n and k exactly.
2. **Estimates distance d.** It uses a randomised information-set search, so d is an *upper bound*.
3. **Screens** on k·d²/n, the standard figure of merit (the gross code scores 12).
4. **Monte Carlo (fast loop).** Codes that pass the screen get thousands of noisy trials decoded with BP-OSD. The model is code capacity: i.i.d. errors with perfect syndromes.
5. **Evolves.** Each torus size keeps its own population. Children come from mutation and crossover, and trivially equivalent codes are skipped.
6. **Circuit level (slow loop, `verify`).** The best survivors get a full noisy syndrome-extraction circuit in Stim, decoded on the detector error model. Each one runs next to the IBM code of the same size as a control.

## Setup (once)

Needs Python 3.10+.

```bash
cd qec-search
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m qec_search.selftest
```

The self-test must say **ALL CHECKS PASSED**. It checks that the five IBM codes come out with their published parameters, that the decoder works, and that the Stim circuits are valid. If anything fails, send Claude the output.

## The loop

```bash
# 1. search (Ctrl-C any time; results so far are kept; re-running the same --out resumes)
python -m qec_search.search --out runs/run1 --hours 2 --seed-known

# 2. circuit-level check of the best finds (+ IBM controls)
python -m qec_search.summarize runs/run1
python -m qec_search.verify runs/run1 --top 5 --ps 0.003

# 3. final report
python -m qec_search.summarize runs/run1
```

Then send Claude **`runs/run1/report_for_claude.md`**, plus `top_candidates.json` if it's asked for.

## Useful knobs

| flag | meaning |
|---|---|
| `--sizes 6x6,12x6,...` | which tori to search (n = 2·l·m) |
| `--terms 3,3` | monomials in A and B. 3,3 gives weight-6 checks like IBM's; try 3,4 or 4,4 |
| `--mixed` | allow mixed monomials x^i·y^j, not just pure powers (less explored) |
| `--workers N` | CPU processes (default: all cores but one) |
| `--ps 0.03,0.05` | noise levels for the fast loop |
| `--shots / --max-fails` | Monte Carlo effort per code |
| `--ler-gate 0.6` | only simulate codes with k·d²/n at least 0.6× the best reference at that size |
| `--dist-trials` | effort for the distance bound (raise for n > 250) |

Rough speed with `ldpc` installed: a few thousand candidates per hour per core. Most candidates have k = 0 and cost almost nothing.

## Files a run produces

- `results.jsonl`: one line per code with k > 0 (genome, n, k, d_ub, Monte Carlo results)
- `zero_k_keys.txt`: codes that turned out to have k = 0, so a resumed run skips them
- `meta.json`: settings, library versions, counters
- `circuit.jsonl`: circuit-level results from `verify`
- `report_for_claude.md`, `top_candidates.json`: written by `summarize`

## Honest limits

- **d is an upper bound.** Finalists need an exact distance computation (e.g. an integer program or a QDistRnd-style search) before anyone claims a parameter.
- **Code-capacity noise is only a screen.** Circuit-level results decide. The circuit uses a generic schedule (all X checks, then all Z checks). That is fair across candidates but worse than IBM's hand-tuned depth-8 schedule, so compare against the controls, not against IBM's published numbers.
- **"Beats reference" means beats 5 IBM codes.** Hundreds of BB and generalised-bicycle codes have been published since 2024, including LLM-guided and RL searches in 2026. Check a find against the literature before calling it new.
- **Equivalence detection covers the obvious symmetries only.** Two different keys can still be the same code.
