# Memory

## Always

## Lessons
- d=16 BB proofs: HiGHS MILP (qec_search.exact_distance) ~25 min; SAT cubing (certify_upstream.py) crawled; use MILP.
- qldpc --open-pr always fails prose pre-flight (draft has unticked equivalence box); write body, commit/push by hand.
- verify/prepush_prose_check.sh needs python3: set QLDPC_PYTHON=<worktree>/.venv/Scripts/python.exe on Windows.
- Local scipy 1.18 HiGHS stalled on 288_12_16 MILP (2% at 36 min) vs 31 min on cloud 1.17; big proofs go to Colab.
- Board shows d<= unless certs/ from upstream verify/certify.py; our MILP is note-level evidence only.
- verify/sat_certify._logicals incomplete-basis bug (our #2273) is fixed upstream by #2343; old certs predate the fix.
- certify_sym 180-18-14 --tlim 600: pending fell to 32 at 70 cpu-h, then rose; len-16 cubes split 17/23 (diverging).
- 288-12-18 bilayer: best GL(2,Z12)+shift+fold start 7.211 (120 checks >7); cold/hot anneal only to 111; not at 7.0.
- Git Bash mangles /mnt/c paths passed to wsl.exe: prefix MSYS_NO_PATHCONV=1 (scripts/run_distqldpc.sh).
- DistQLDPC one-sided (export_distqldpc.py: Hz:=I, Gx:=0 row) is 4.6x faster than its symplectic form on 144-12-12.
- Exact distances go through qec_search/distqldpc.solve_many only; a script-local DistQLDPC call lost the orbital split.
- Never `uv run --python 3.11`: its managed CPython is blocked from executing here (os error 5) and new venvs adopt it.
- Hunt screens use qec_search/frontier.screen_distance (stop_at = bar - 1); a script-local distance bound made screens ~10x slower.
- Never run two submit.py at once: both `git fetch upstream` and the ref update collides ("incorrect old value").
- A wait loop testing `ps | grep "[x]yz"` matches its own `bash -c` line if that line contains xyz; poll a log file instead.
- graft_reduce.py: a 2k+20k RIS filter passed a d-losing graft (392->375 read d<=16 at 200k); keep 20k + 2x200k.
- gf2_fast RIS: small consecutive seeds (1,2,3..) miss rare light logicals; hash seeds (zlib.crc32) as ris_ladder does.
- Windows python print/write_text emits CRLF; a .sh or command file fed to wsl.exe then carries \r. Write bytes with \n.
- wsl.exe inside a 'while read' loop eats the loop's stdin; always run it with < /dev/null.
