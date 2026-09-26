"""Submit one code to the qLDPC Challenge from this fork.

    python lab/scripts/submit.py bb_216_8_16_18x6 --handle @yourhandle            # dry run (default)
    python lab/scripts/submit.py bb_216_8_16_18x6 --handle @yourhandle --for-real # branch, commit, push, open PR

The lab lives on our fork's main, so a submission must never branch from main. Instead:
  1. git fetch upstream
  2. a worktree at <fork>/.worktrees/<name>, detached at upstream/main (clean upstream tree)
  3. inside it, the challenge's own `qldpc submit` (verifier + prose gate), with --dry-run or --open-pr.
     --open-pr creates branch submit-<n>-<k>-<d>, commits only codes/ (+ note), pushes to origin
     (our fork) and opens the PR against unitaryfoundation/qldpc-challenge.

Refuses candidates whose manifest status is "do_not_submit" or "hold" unless --force.
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
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[1]
FORK = LAB.parent
UPSTREAM_REPO = "unitaryfoundation/qldpc-challenge"
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


def qldpc_args(info: dict, npz: Path, handle: str, *, family: str, model: str, has_layout: bool,
               circuits: bool, note: Path | None, for_real: bool) -> list[str]:
    """Arguments for the challenge CLI's `submit` subcommand."""
    cmd = ["submit", str(npz), "--authors", handle, "--family", family, "--construction", info["construction"]]
    if model:
        cmd += ["--model", model]
    if has_layout:
        cmd += ["--layers", str(info.get("layers", 2)), "--coords", str(npz)]
    if not circuits:
        cmd.append("--no-circuit")
    if note:
        cmd += ["--note-file", str(note)]
    cmd.append("--open-pr" if for_real else "--dry-run")
    return cmd


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
    args = ap.parse_args(argv)

    manifest_path = Path(args.manifest)
    with open(manifest_path) as f:
        info = json.load(f)[args.name]
    if reason := refusal(args.name, info, args.force):
        sys.exit(reason)

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
    cmd = qldpc_args(info, npz, args.handle, family=args.family, model=args.model, has_layout=has_layout,
                     circuits=args.circuits, note=note, for_real=args.for_real)
    r = sh([uv, "run", "--quiet", "python", "cli/qldpc.py", *cmd], cwd=wt, check=False, env=env)
    if not args.for_real:
        print(f"\nDry run only (checked against upstream/main in {wt}). Re-run with --for-real to open the PR.")
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
