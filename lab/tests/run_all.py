"""Run the test suite without pytest:  python tests/run_all.py"""
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_core  # noqa: E402

failed = 0
for name in sorted(n for n in dir(test_core) if n.startswith("test_")):
    t = time.time()
    try:
        getattr(test_core, name)()
        print(f"PASS  {name}  ({time.time() - t:.1f}s)")
    except Exception:
        failed += 1
        print(f"FAIL  {name}")
        traceback.print_exc()
print("\nall tests passed" if not failed else f"\n{failed} test(s) failed")
sys.exit(1 if failed else 0)
