# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Submit one code to the qLDPC Challenge from this fork.

    python lab/scripts/submit.py bb_216_8_16_18x6 --handle @yourhandle            # dry run (default)
    python lab/scripts/submit.py bb_216_8_16_18x6 --handle @yourhandle --for-real # branch, commit, push, open PR

The lab lives on our fork's main, so a submission must never branch from main. Instead:
  1. git fetch upstream
  2. a worktree at <fork>/.worktrees/<name>, detached at upstream/main (clean upstream tree)
  3. inside it, the challenge's own `qldpc submit --dry-run --json` builds the doc, and upstream's
     verify/validate_candidate.py gates it (verify + refute + board dedup); a failed gate stops --for-real.
  4. then `qldpc submit` again (verifier + prose gate), with --dry-run or --open-pr.
     --open-pr creates branch submit-<n>-<k>-<d>, commits only codes/ (+ note), pushes to origin
     (our fork) and opens the PR against unitaryfoundation/qldpc-challenge.

Refuses candidates whose manifest status is "do_not_submit" or "hold" unless --force, and always
refuses one without a finished exact-distance proof of its genome in results/proofs/, or whose
witnessed d (what the board will show) differs from the proven d.
Runs the challenge CLI through `uv` (installed in the lab venv), so no bash is needed.
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[1]
FORK = LAB.parent
UPSTREAM_REPO = "unitaryfoundation/qldpc-challenge"
PROOFS = LAB / "results/proofs"
BLOCKED = ("hold", "do_not_submit")


def sh(cmd, cwd, check=True, env=None):
    print(f"$ {shlex.join(cmd)}", flush=True)
    return subprocess.run(cmd, cwd=cwd, check=check, env=env)


def refusal(name: str, info: dict, force: bool) -> str | None:
    """Reason to refuse this candidate, or None if it may be submitted."""
    status = info.get("status", "submit")
    if status in BLOCKED and not force:
        return f"refusing: {name} has status '{status}': {info.get('status_reason')}  (use --force to override)"
    return None


def params(info: dict) -> tuple[int, int, int]:
    """(n, k, d) from a manifest entry's "[[n,k,d]]"."""
    n, k, d = (int(x) for x in info["params"].strip("[]").split(","))
    return n, k, d


def _genome_key(g: dict) -> tuple:
    """A BB genome (l, m, A, B) or a tile code (qec_search.tile: box B, bulk l x m, X-tile edges)."""
    if "tile" in g:
        return "tile", g["B"], g["l"], g["m"], tuple(sorted(map(tuple, g["tile"])))
    return g["l"], g["m"], tuple(sorted(map(tuple, g["A"]))), tuple(sorted(map(tuple, g["B"])))


def proof_refusal(name: str, info: dict, proofs_dir: Path) -> str | None:
    """Reason this candidate lacks a finished exact-distance proof of its own genome and d, or None."""
    n, k, d = params(info)
    path = proofs_dir / f"distance_proof_{n}_{k}_{d}.json"
    if not path.exists():
        return f"no distance proof for {name}: expected {path}"
    proof = json.loads(path.read_text(encoding="utf-8"))
    if "proven exact" not in proof.get("result", ""):
        return f"{path.name} is not a finished proof: {proof.get('result')}"
    if proof.get("d") != d:
        return f"{path.name} proves d={proof.get('d')}, the manifest claims d={d}"
    if _genome_key(proof["genome"]) != _genome_key(info["genome"]):
        return f"{path.name} is for a different genome than {name}"
    return None


def ladder_refusal(info: dict, log: Path, min_trials: int = 5_000_000, min_rungs: int = 3) -> str | None:
    """Reason the candidate's RIS ladder (scripts/ris_ladder.py) does not back its claimed d as an upper bound, or
    None. The manifest's "ladder_npz" names the code's rows in the log; it needs at least min_rungs fresh-seed rungs,
    one of at least min_trials, and the lightest logical any rung found must be exactly the claimed d. This backs a
    d<= board claim only (user decision 15 Oct 2026); it is not a proof."""
    _, _, d = params(info)
    key = info.get("ladder_npz")
    if not key:
        return "manifest entry has no ladder_npz"
    if not log.exists():
        return f"no ladder log {log}"
    rungs = [r for r in (json.loads(x) for x in log.read_text(encoding="utf-8").splitlines() if x.strip())
             if r["npz"] == key]
    if len({r["seed"] for r in rungs}) < min_rungs:
        return f"{key}: {len(rungs)} rung(s) in {log.name}, need {min_rungs} with distinct seeds"
    if max(r["trials"] for r in rungs) < min_trials:
        return f"{key}: largest rung {max(r['trials'] for r in rungs)} trials, need {min_trials}"
    if (low := min(r["d_ub"] for r in rungs)) != d:
        return f"{key}: the ladder's lightest logical has weight {low}, the manifest claims d={d}"
    return None


def wsl_path(p: str | Path) -> str:
    """C:\\a\\b -> /mnt/c/a/b; anything that is not an absolute Windows path is returned unchanged."""
    s = str(p)
    if len(s) > 2 and s[1] == ":" and s[2] in "\\/":
        return f"/mnt/{s[0].lower()}/" + s[3:].replace("\\", "/")
    return s


def wsl_uv_run(wt: Path, args: list[str]) -> list[str]:
    """`uv run python <args>` in the worktree under WSL, where gf2_fast is built, so qldpc's witness search reaches
    the distances the ladder measured. The challenge needs Python >= 3.12 (WSL's system python is 3.10), so uv runs
    a managed 3.12 with gf2_fast built for it once:
        uv run --python 3.12 --with pybind11 --with setuptools python verify/setup_gf2_fast.py build_ext             --build-lib ~/build/gf2_fast312 --build-temp /tmp/gf2tmp312
    The WSL env lives outside the worktree so it never collides with the Windows .venv."""
    inner = ("cd {} && UV_PROJECT_ENVIRONMENT=$HOME/.venvs/qldpc-wsl PYTHONPATH=$HOME/build/gf2_fast312 PYTHONUTF8=1 "
             "$HOME/.local/bin/uv run --python 3.12 --quiet python {}").format(shlex.quote(wsl_path(wt)),
                                                                 shlex.join(wsl_path(a) for a in args))
    return ["wsl.exe", "-e", "bash", "-lc", inner]


def claim_mismatch(doc: dict, d: int) -> str | None:
    """Reason the distance qldpc would submit differs from the proven one, or None."""
    got = doc.get("distance", {}).get("d")
    return None if got == d else f"qldpc would submit d={got}, but the proven distance is {d}"


def free_slug(codes_dir: Path, n: int, k: int, d: int) -> str:
    """First unused board slug for these parameters: n-k-d, then n-k-d-b, -c, ... (the board's own convention)."""
    base = f"{n}-{k}-{d}"
    for suffix in [""] + [f"-{c}" for c in "bcdefghijklmnopqrstuvwxyz"]:
        if not (codes_dir / f"{base}{suffix}.json").exists():
            return base + suffix
    raise RuntimeError(f"no free slug for {base}")


def qldpc_args(info: dict, npz: Path, handle: str, *, family: str, model: str, has_layout: bool,
               circuits: bool, note: Path | None, for_real: bool, json_doc: bool = False) -> list[str]:
    """Arguments for the challenge CLI's `submit` subcommand."""
    cmd = ["submit", str(npz), "--authors", handle, "--family", family, "--construction", info["construction"]]
    if model:
        cmd += ["--model", model]
    if info.get("provenance_notes"):
        cmd += ["--notes", info["provenance_notes"]]
    if has_layout:
        cmd += ["--layers", str(info.get("layers", 2)), "--coords", str(npz)]
    if not circuits:
        cmd.append("--no-circuit")
    if note:
        cmd += ["--note-file", str(note)]
    cmd.append("--open-pr" if for_real else "--dry-run")
    if json_doc and not for_real:
        cmd.append("--json")
    return cmd


def doc_from_dry_run(stdout: str) -> dict:
    """The submission doc from `qldpc submit --dry-run --json`: the current CLI prints one JSON object on stdout with
    the doc under "doc"; older ones printed it after a "would write" line."""
    try:
        res = json.loads(stdout)
    except json.JSONDecodeError:
        res = None
    if isinstance(res, dict):
        if isinstance(res.get("doc"), dict):
            return res["doc"]
        raise ValueError(f"no submission doc in the dry-run result: {res.get('error') or res.get('stage')}")
    marker = stdout.find("--dry-run: would write")
    start = stdout.find("\n{", marker)
    if marker < 0 or start < 0:
        raise ValueError("no submission JSON in the dry-run output")
    return json.loads(stdout[start + 1:])


def gate_verdict(verdict: dict) -> str | None:
    """Reason validate_candidate's verdict blocks a submission, or None if it passed."""
    if verdict.get("passed"):
        return None
    return "validate_candidate failed: " + ("; ".join(verdict.get("labels", [])) or "no labels")


def prepare_worktree(name: str, for_real: bool) -> Path:
    """A worktree detached at upstream/main; refreshed if unused, refused if it already holds a submission."""
    sh(["git", "fetch", "upstream", "--quiet"], cwd=FORK)
    wt = FORK / ".worktrees" / name
    if not wt.exists():
        sh(["git", "worktree", "add", "--quiet", "--detach", str(wt), "upstream/main"], cwd=FORK)
        return wt
    on_branch = subprocess.run(["git", "symbolic-ref", "-q", "HEAD"], cwd=wt, capture_output=True).returncode == 0
    if on_branch:
        if for_real:
            sys.exit(f"{wt} already holds a submission branch. Amend it there, or remove it with "
                     f"`git worktree remove {wt}` once the PR is merged or closed.")
        return wt
    sh(["git", "checkout", "--quiet", "--detach", "upstream/main"], cwd=wt)
    return wt


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", help="candidate name, e.g. bb_216_8_16_18x6")
    ap.add_argument("--handle", required=True, help="your GitHub handle, e.g. @fullymiddleaged")
    ap.add_argument("--manifest", default=str(LAB / "results/candidates/manifest.json"))
    ap.add_argument("--model", default="Claude Opus 5.5", help="exact AI model version credited on the entry ('' to omit)")
    ap.add_argument("--family", default="bivariate-bicycle")
    ap.add_argument("--note", type=Path, help="research note (notes/<n>-<k>-<d>.md) to ship with the code")
    ap.add_argument("--circuits", action="store_true", help="let the challenge tool run its circuit simulation (slow)")
    ap.add_argument("--for-real", action="store_true", help="branch, commit, push and open the PR (default: dry run)")
    ap.add_argument("--force", action="store_true", help="ignore a hold / do_not_submit status")
    ap.add_argument("--upper-bound", type=Path, metavar="LADDER_JSONL",
                    help="submit a d<= claim backed by this RIS ladder log instead of an exact proof")
    ap.add_argument("--wsl", action="store_true",
                    help="run the challenge CLI and gate in WSL, where gf2_fast is built (needed at large n)")
    args = ap.parse_args(argv)

    manifest_path = Path(args.manifest)
    with open(manifest_path) as f:
        info = json.load(f)[args.name]
    if reason := refusal(args.name, info, args.force):
        sys.exit(reason)
    # our aim is proven distances: no finished proof of this genome at this d, no submission, unless the user
    # chose an upper-bound claim for this code and a full fresh-seed ladder backs it
    if args.upper_bound:
        reason = ladder_refusal(info, args.upper_bound.resolve())
    else:
        reason = proof_refusal(args.name, info, PROOFS)
    if reason:
        sys.exit(f"refusing: {reason}")

    npz = (manifest_path.parent / f"{args.name}.npz").resolve()
    with np.load(npz) as z:
        has_layout = "coords" in z.files
    note = args.note.resolve() if args.note else None
    if args.for_real and not note:
        print("warning: no --note; CONTRIBUTING.md asks for a research note with each submission", flush=True)

    # uv lives in the lab venv next to this interpreter; the challenge CLI runs in its own uv env
    # (PYTHONUTF8: the challenge CLI prints "≤", which a Windows cp1252 console can't encode)
    env = dict(os.environ, PYTHONUTF8="1",
               PATH=os.pathsep.join([str(Path(sys.executable).parent), os.environ.get("PATH", "")]))
    uv = shutil.which("uv", path=env["PATH"])
    if not uv:
        sys.exit("uv not found: pip install uv in the lab venv")
    if args.for_real:
        # gh pr create in a fork needs to know the base repo; this is local git config
        sh(["gh", "repo", "set-default", UPSTREAM_REPO], cwd=FORK)

    wt = prepare_worktree(args.name, args.for_real)
    # qldpc submit always writes codes/<n>-<k>-<d>.json; when the board already has that slug (another code with
    # these parameters, e.g. IBM's gross code for our bilayer layout of it) it must go in as <slug>-b by hand
    if (slug := free_slug(wt / "codes", *params(info))) != "{}-{}-{}".format(*params(info)):
        sys.exit(f"codes/{'{}-{}-{}'.format(*params(info))}.json already exists upstream. Run `qldpc submit ... --out <tmp>` "
                 f"in {wt}, then commit it as codes/{slug}.json with notes/{slug}.md (as PR #2394 did).")
    common = dict(family=args.family, model=args.model, has_layout=has_layout, note=note)

    # gate: upstream's trusted validator on the exact doc qldpc would write (no circuits: they don't affect it)
    def uv_run(script_args: list[str]) -> list[str]:
        return wsl_uv_run(wt, script_args) if args.wsl else [uv, "run", "--quiet", "python", *script_args]

    cmd = qldpc_args(info, npz, args.handle, **common, circuits=False, for_real=False, json_doc=True)
    print(f"$ uv run python cli/qldpc.py {shlex.join(cmd)}" + ("  (in WSL)" if args.wsl else ""), flush=True)
    dry = subprocess.run(uv_run(["cli/qldpc.py", *cmd]), cwd=wt, env=env,
                         capture_output=True, text=True, encoding="utf-8")
    if dry.returncode:
        print(dry.stdout, dry.stderr, sep="\n")
        sys.exit(dry.returncode)
    doc = doc_from_dry_run(dry.stdout)
    if reason := claim_mismatch(doc, params(info)[2]):
        sys.exit(f"refusing: {reason}")
    with tempfile.TemporaryDirectory() as tmp:
        doc_path = Path(tmp) / f"{args.name}.json"
        doc_path.write_text(json.dumps(doc), encoding="utf-8")
        print(f"$ uv run python verify/validate_candidate.py {doc_path}", flush=True)
        val = subprocess.run(uv_run(["verify/validate_candidate.py", str(doc_path)]),
                             cwd=wt, env=env, capture_output=True, text=True, encoding="utf-8")
    try:
        verdict = json.loads(val.stdout)
    except json.JSONDecodeError:
        print(val.stdout, val.stderr, sep="\n")
        sys.exit("validate_candidate produced no verdict")
    print(f"validate_candidate: passed={verdict['passed']}  " + " | ".join(verdict.get("labels", [])), flush=True)
    if reason := gate_verdict(verdict):
        if args.for_real:
            sys.exit(reason)
        print(f"warning: {reason} (a --for-real run would stop here)", flush=True)

    cmd = qldpc_args(info, npz, args.handle, **common, circuits=args.circuits, for_real=args.for_real)
    if args.wsl and args.for_real:
        sys.exit("--wsl stops after the gate: write the file with `qldpc submit ... --out codes` under WSL in "
                 f"{wt}, then commit, push and open the PR by hand (git and gh are set up on Windows only).")
    r = sh(uv_run(["cli/qldpc.py", *cmd]), cwd=wt, check=False, env=env)
    if not args.for_real:
        print(f"\nDry run only (checked against upstream/main in {wt}). Re-run with --for-real to open the PR.")
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
