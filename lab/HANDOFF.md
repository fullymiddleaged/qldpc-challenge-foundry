# Handoff: qldpc-lab

> **Layout changed on 2026-09-26.** The lab now lives in `lab/` on the `main` branch of the fork
> (fullymiddleaged/qldpc-challenge-foundry), not in a separate repo next to it. §1 and §2 describe the
> original two-repo plan. `lab/CLAUDE.md`, `lab/scripts/setup_workspace.sh` and `lab/scripts/submit.py`
> describe the current setup. Colab notebooks are in `lab/colab/`.

*Written 2026-09-26 at the end of the first working sessions (Claude in the Claude app, runs on Google Colab). This is for picking the work up in VS Code with Claude Code.*

---

## 0. TL;DR

- We built a **search → layout → certify → submit** pipeline for quantum LDPC codes of the bivariate bicycle (BB) family: the family behind IBM's [[144,12,12]] "gross code".
- The live frontier we compete on is the **Unitary Foundation qLDPC Challenge** leaderboard (github.com/unitaryfoundation/qldpc-challenge): an auto-verified, Pareto-ranked board with cells by locality (2D single-layer / 2D bilayer / unrestricted) × check weight.
- Our pure code search found nothing new in the crowded *unrestricted* cells. But once our codes have a **2D bilayer layout** they beat the current entries in the **2D-local bilayer, weight ≤ 6** cell:

| Code | Layout radius (cap 7.0) | Distance | Status |
|---|---|---|---|
| [[216,8,16]] | 7.000 | **proven exact** (MILP, 23 min) | ready to submit |
| [[288,12,16]] | 7.000 | **proven exact** (MILP, 31 min) | ready to submit |
| [[288,8,20]] | 6.708 | upper bound 20 only | **hold**: exact proof pending |
| [[144,8,12]] | 6.708 | upper bound | optional (weak: displaced by anyone laying out the gross code) |
| IBM gross [[144,12,12]] | 6.708 | known exact 12 | **do not submit** (IBM's code) |

- The challenge's own checker (dry run, 2026-09-26) verified all five: OK, locality `local-2d-bilayer`, witness distances matching ours. Log: `results/logs/checker_bilayer_candidates_dryrun_2026-09-26.txt`.
- **Honest framing:** these are records *under the board's rules*. The board's "bilayer" metric ignores ancilla qubits and wiring, so it is much easier than real chip layout (a 2026 npj QI paper needs 5 tiers of couplers for the gross code). Real hardware-grade layout is the obvious research step up (Milestone 4 below).

---

## 1. Setup after you fork

### 1.1 Fork and unpack
1. On GitHub, fork `unitaryfoundation/qldpc-challenge` to your account.
2. Make a workspace folder and unzip this handoff into it, so you get:
   ```
   qldpc-workspace/
     qldpc-lab/          <- this repo (from the zip)
     qldpc-challenge/    <- created by the setup script (your fork)
   ```
3. From `qldpc-workspace/qldpc-lab`, run:
   ```bash
   bash scripts/setup_workspace.sh <your-github-username>
   ```
   This does six things. It clones your fork next to this repo and adds the `upstream` remote. It creates `.venv` and installs this package with `stim`, `ldpc` and `pytest`. It installs `uv`, which the challenge tool needs. It runs `python -m qec_search.selftest`, which checks the IBM codes, the decoder and the Stim circuits. It runs `pytest`. Finally it smoke-tests `./qldpc targets` in the fork, and does `git init` for qldpc-lab if needed.
4. Create an empty GitHub repo called `qldpc-lab` (private is fine to start). Then run `git remote add origin … && git push -u origin main`.
5. Open the workspace in VS Code with `code qldpc.code-workspace`. Both repos show in the sidebar.

**Windows:** run the script from **Git Bash** or **WSL**; the challenge's `./qldpc` launcher is a bash script. In plain PowerShell, do the steps by hand: `python -m venv .venv; .venv\Scripts\activate; pip install -e ".[sim,dev]"; pip install uv`.

### 1.2 First prompt for Claude Code
Paste this once the workspace is open:

> Read CLAUDE.md and HANDOFF.md. Run the test suite and `python -m qec_search.selftest`, and report anything that fails. Then check the two "ready to submit" codes with `python scripts/submit.py bb_216_8_16_18x6 --handle @<me>` (dry run) and the same for `bb_288_12_16_24x6`. Also show me the current bilayer weight-≤6 frontier from `../qldpc-challenge/qldpc targets`. Don't submit anything until I say so.

---

## 2. How the two repos fit together

- **qldpc-lab (ours)** holds all code, experiments, results, proofs and docs. It is the long-lived project, and it is where "advancing the field" happens.
- **qldpc-challenge (your fork)** is used only for its verifier (`./qldpc`) and for submission branches. **Never commit our code into the fork's `main`.** The challenge wants one JSON per pull request.
- **Submission flow:** `scripts/submit.py <candidate> --handle @you --for-real` does the following:
  1. fetches `upstream`;
  2. creates the branch `submit/<candidate>` from `upstream/main`;
  3. runs `./qldpc submit` with our construction note, layout flags and model credit;
  4. commits `codes/<n>-<k>-<d>.json`, pushes it to your fork, and prints the "open PR" link (or uses `gh` with `--gh`).

  I tested this end to end against a stand-in local fork (branch, commit, push, link). I have **not** run it against real GitHub.
- **Contributing tools upstream** (e.g. our MILP distance certifier, see Milestone 2) is a separate kind of PR. Do it on a feature branch of the fork cut from `upstream/main`, after discussing it with the maintainers in an issue.

---

## 3. What is in the package

| Module | What it does | Tested? |
|---|---|---|
| `gf2.py` | GF(2) rank, RREF, nullspace, solve | yes |
| `bbcode.py` | BB code from (l, m, A, B). `distance_upper_bound` (randomised information sets, *upper bound*). `components()` detects split codes. `signature()` is an eigenvalue fingerprint for relabelled duplicates. `canonical_key()` covers the simple symmetries only. | yes |
| `known_codes.py` | IBM's 5 reference codes (Bravyi et al., Nature 2024) plus their signatures | yes |
| `decoders.py` | BP-OSD via `ldpc` (v2 API), numpy BP + OSD-0 fallback | yes (ldpc on Colab) |
| `evaluate.py` | Code-capacity Monte Carlo (fast loop). `circuit_level_ler` (slow loop, Stim) | fast loop yes, circuit loop only via the self-test |
| `circuits.py` | Stim memory-Z circuit for any CSS code (generic schedule: all X checks, then all Z), plus DEM → matrices | self-test stage 3 checks determinism and the DEM. The Colab self-test ran, but its stage-3 output was never reviewed. **Circuit-level LERs have never been produced or compared.** |
| `search.py` | Evolutionary search with per-torus-size niches, resumable JSONL output, dedupe by key + signature, split-code filter | yes (Colab: ~19k candidates/h on 2 cores) |
| `summarize.py` | Compact markdown report + `top_candidates.json` | yes |
| `verify.py` | Circuit-level check of top candidates, with IBM controls | untested with real Stim |
| `layout.py` | Bilayer layout: folded-torus start + simulated annealing on a grid, cap 2 per site, minimising the max check diameter | yes |
| `export.py` | Pareto-front export to `.npz` (hx, hz[, coords]) + manifest; skips IBM codes; `--layout` | yes |
| `exact_distance.py` | **MILP proof of minimum distance** (HiGHS via scipy). Uses translation symmetry and d_X = d_Z. Reports a proven lower bound on timeout. | yes: reproduces 6, 10, 12 for the IBM codes |
| `selftest.py` | 3-stage self-test | yes |

`notebooks/qec_search_colab.ipynb` is the Colab version, with everything embedded. Use it for multi-hour jobs.

---

## 4. Findings so far

### 4.1 Search (baseline, Colab, 2026-09-25)
- The run checked 19,234 IBM-style candidates (weight 6, pure-power terms) in 1 h. 3,618 had k > 0, and 225 got the Monte Carlo test.
- It rediscovered IBM's [[72,12,6]], [[90,8,10]], [[108,8,10]] and [[144,12,12]] unaided, and nothing beat them on k·d²/n.
- The IBM-style space is small. There are only **1,285 / 3,546 / 16,575 / 2,298** distinct codes on the 6×6 / 9×6 / 12×6 / 15×3 tori, so it is exhaustible and almost certainly already mined. The `--mixed` space (x^i·y^j terms) is millions of times larger. **The `mixed1` run has not been done yet.**
- **Traps we hit:**
  - The "best" [[144,12,12]] was the gross code relabelled; graph isomorphism confirmed it.
  - [[288,24,12]] was two gross codes side by side.
  - About 55% of k > 0 codes in a test run were copies of this kind. The signature and components filters now handle them.

### 4.2 Leaderboard (snapshot 2026-09-26, from `./qldpc targets --n 300`)
- There are 1,443 codes across 12 cells, and the board is active daily.
- **Unrestricted, w ≤ 6, n ≤ 300, top by k·d²/n:** [[288,16,16]] 14.22, [[254,14,16]] 14.11, [[288,12,18]] 13.50, [[248,10,18]] 13.06, [[192,12,14]] 12.25, [[210,10,16]] 12.19, [[252,12,16]] 12.19, [[144,12,12]] 12.00, [[192,16,12]] 12.00, [[270,8,20]] 11.85. All of our codes are dominated here.
- **2D-local bilayer, w ≤ 6, n ≤ 300:** [[300,8,20]] 10.67, [[280,6,22]] 10.37, [[264,8,18]] 9.82, [[252,12,14]] 9.33, [[168,6,16]] 9.14, [[192,12,12]] 9.00, [[240,8,16]] 8.53, [[238,6,18]] 8.17, [[154,6,14]] 7.64, [[108,8,10]] 7.41.
  - Our [[288,8,20]] dominates [[300,8,20]].
  - Nothing there dominates [[288,12,16]], [[216,8,16]] or [[144,8,12]]. Any code that dominates ours must also score higher on k·d²/n, so it would have appeared in this list.
- **Unrestricted, w ≤ 8:** the top is 38.48 ([[232,62,12]]), and [[168,20,14]] is 23.33. Weight-8 codes are far stronger.
- Full logs are in `results/logs/`.

### 4.3 Board rules we verified (TRACKS.md / SCHEMA.md / cli/qldpc.py)
- **Locality classes:** `local-2d-single` (radius ≤ 4.0), `local-2d-bilayer` (≤ 7.0, up to 2 layers) and `unrestricted` (no layout). The classes nest.
- **Interaction radius** = the claimed maximum *check diameter*: the largest distance between any two data qubits of one check. Ancillas are not modelled.
- **Placement:** at most `layers` qubits per site, and distinct sites ≥ 1 apart. n ≤ 700, and ≤ 600 coordinates.
- **Input:** `.npz` with `hx`/`hz` (also `HX`, `H_X`, …), plus optional `coords`. Flags: `--layers`, `--coords`, `--family`, `--construction`, `--model`, `--no-circuit`/`--circuits`, `--dry-run`, `--open-pr`.
- **Distance:** a witness logical proves an **upper bound**. "Exact" needs the maintainers' server certification. Non-records are accepted, with the record status shown for information.

### 4.4 Literature context
- *Placing and routing quantum LDPC codes in multilayer superconducting hardware* (npj QI 2026, arXiv 2507.23011) needs a nearest-neighbour tier plus **4 more tiers** to keep the gross code's couplers short. It routes real couplers, including ancillas, and it is the right comparison for Milestone 4.
- Recent code discovery with LLMs/RL is active: IBM's LLM-evolved BB codes (2026), OmniQEC (arXiv 2607.25865), and multi-agent qLDPC discovery. Check these before calling a code "new".

---

## 5. Open threads (do these first)

1. **Submit [[216,8,16]] and [[288,12,16]]**, which are proven and dry-run verified: `python scripts/submit.py <name> --handle @you --for-real`, one at a time. Watch the CI checks and reply to the maintainers.
2. **Prove [[288,8,20]].** Run the MILP for 10+ h on a machine that won't restart (Colab cell 7, or locally with `nohup`). Outcomes:
   - "PROVEN exactly 20": set the manifest status to `submit`, add the proof sentence, and submit.
   - Only "≥ X": your call. It can go in labelled as an upper bound, which is normal on the board.
   - The solver finds a logical of weight < 20: the code is weaker than claimed; do not submit it as 20.

   The first attempt died after 2h20m in a host restart. Splitting the problem per logical is not faster.
3. **Run the `mixed1` search** (`--mixed`, 3 h+), then `export --layout` and dry-run the Pareto front.
4. **Validate the circuit-level path.** Run `python -m qec_search.selftest` and confirm stage 3 prints OK for the noiseless-circuit and detector-error-model lines. That output has never been reviewed, and no circuit-level error rates have been produced yet. Run `verify` on the gross code with low shots first, and compare the order of magnitude with published numbers.

---

## 6. Roadmap: where the research frontier meets the leaderboard

Every milestone should produce something that is both **useful to the field** (a tool, a method, or a clearly-stated result) and **measurable on the board**.

### M1. Layout-aware code search (co-design)
- **Idea:** stop searching codes first and laying them out afterwards. Put layout feasibility inside the search: restrict the monomial offsets so that a folded-torus embedding keeps check diameters ≤ 7 (or ≤ 4 for single-layer), and score candidates on (k, d) given that constraint.
- **Why:** the bilayer frontier at n ≤ 300 tops out at k·d²/n ≈ 11, while the unrestricted frontier is at 14+. Codes that are both strong and local are exactly what hardware needs, and the board rewards them directly.
- **Deliverables:** a `--local {bilayer,single}` search mode; a fast analytic layout bound for pre-screening; new records in bilayer w ≤ 6 and w ≤ 8, and an attempt at the much harder single-layer cells.
- **Done when:** several new non-dominated bilayer entries, each with an exact-distance proof.

### M2. Certification at scale
- **Idea:** the board's distances are mostly upper bounds. Make exact certification cheap. Use stronger symmetry breaking (the full automorphism group, not just translations), and add lower-bound methods (MILP dual bounds; Brouwer–Zimmermann-style for BB codes).
- **Why:** it separates solid entries from lucky ones, and it can be contributed upstream as the "server certification" step the challenge describes.
- **Deliverables:** `certify` CLI and batch mode; a report of which frontier entries are exact vs upper bounds; an upstream issue/PR offering the tool.
- **Done when:** [[288,8,20]] is certified in under 2 h, and the method is benchmarked on 10+ board entries.

### M3. Circuit-level truth
- **Idea:** k·d²/n is a proxy. Validate `circuits.py`, then produce circuit-level logical error rates per round for our records, with IBM controls. Optimise the syndrome-extraction schedule (coloration circuits), not just the generic X-then-Z one.
- **Why:** this is what hardware groups actually compare; the board's tool also has a `--circuits` option.
- **Done when:** our 288-qubit records have LER-vs-p curves next to the gross code under identical noise.

### M4. Hardware-realistic layout (research-grade)
- **Idea:** add ancillas, physical couplers, crossings and tiers to the layout model, like the npj QI 2026 paper. Optimise layouts (SA / ILP / force-directed) under realistic constraints.
- **Benchmark:** the gross code needs 1 + 4 tiers there. Matching that with fewer tiers, or shorter couplers, would be a publishable result.
- **Also:** propose to the challenge a stricter locality track that includes ancillas.

### M5. Beyond BB codes
- **Families:** generalized bicycle (GB), two-block group algebra (2BGA), coprime-BB, and weight-8 codes (the unrestricted w ≤ 8 frontier is far higher). Reuse the same pipeline: construction → dedupe → layout → certify → submit.

### M6. Decoders (longer term)
- **Idea:** real-time qLDPC decoding (BP-OSD speed/accuracy, Relay-BP, deadline-aware decoding) was the original "hot problem". With circuit-level infrastructure from M3, this becomes approachable.

---

## 7. Troubleshooting

| Problem | Fix |
|---|---|
| `pip install stim ldpc` fails | Needs Python 3.10–3.13 with wheels. Use `python3.12 -m venv .venv`. |
| `./qldpc: Permission denied` / not found | `chmod +x ../qldpc-challenge/qldpc`. On Windows use Git Bash/WSL. |
| `./qldpc` complains about uv | `pip install uv` inside the venv, then retry. |
| Self-test stage 3 is slow | It decodes 2,000 circuit samples; a few minutes is normal. |
| MILP never finishes | Expected for d ≥ 20 at n = 288. Run it in the background with `--hours`, and use the proven lower bound. |
| `submit.py` refuses a code | The manifest status is `hold` / `do_not_submit`. Fix the reason; use `--force` only deliberately. |
| Layout radius stuck just above 7 | Rerun `find_layout` with `pad=5` and more steps/restarts; the pad gives qubits room. |

---

## 8. Sources
- qLDPC Challenge: https://github.com/unitaryfoundation/qldpc-challenge (TRACKS.md, schema/SCHEMA.md, cli/qldpc.py); leaderboard https://unitaryfoundation.github.io/qldpc-challenge/
- Bravyi et al., *High-threshold and low-overhead fault-tolerant quantum memory*, Nature 627, 778 (2024): the IBM BB codes
- *Placing and routing quantum LDPC codes in multilayer superconducting hardware*, npj Quantum Information (2026), arXiv 2507.23011
- *Tour de gross: a modular quantum computer based on bivariate bicycle codes* (IBM)
- IBM Research blog, *Can LLMs discover quantum error correction codes?* (2026); OmniQEC, arXiv 2607.25865
- Error Correction Zoo, bivariate bicycle codes: https://errorcorrectionzoo.org/c/qcga
