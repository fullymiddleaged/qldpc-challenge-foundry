# Colab notebooks

Heavy jobs go here: multi-hour searches, long MILP distance proofs, circuit-level runs. Colab gives you
more cores than a laptop, and a runtime that survives closing your editor.

## How the notebooks get the code
Each notebook's first cell clones **this fork** (`fullymiddleaged/qldpc-challenge-foundry`) at the branch you
pick (`BRANCH`, default `main`) and runs `pip install -e lab[sim]`. The notebook therefore runs whatever is
**pushed** to that branch: push lab changes before you start a Colab run. The cell prints the commit it
is using; copy that into any result you bring back.

Don't embed code in a notebook (the old handoff notebook carried a base64 zip). It drifts from `lab/qec_search`.

## Getting results back
Runs write to Google Drive (`MyDrive/qec-runs/<run>`) so they survive a runtime reset. To keep a result:
download it (the last cell does this), put it under `lab/results/` (`proofs/`, `logs/`, `candidates/`), and commit it.
Raw run directories stay out of git (`runs/` is ignored).

## Notebooks
| Notebook | Use it for |
|---|---|
| `qec_search_colab.ipynb` | search → summarise → circuit check → export + layout → dry-run against the challenge checker → MILP distance proof |

Open one in Colab from GitHub: `https://colab.research.google.com/github/fullymiddleaged/qldpc-challenge-foundry/blob/main/lab/colab/<notebook>.ipynb`
