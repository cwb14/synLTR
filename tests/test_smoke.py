import os
import subprocess
import sys


def test_entry_sets_thp_guard():
    """The console entry sets NUMPY_MADVISE_HUGEPAGE=0 before numpy loads unless the user chose a value, which then
    stands; checked in a fresh interpreter without the variable (the test session's conftest sets it first)."""
    env = {k: v for k, v in os.environ.items() if k != "NUMPY_MADVISE_HUGEPAGE"}
    code = ("import os, sys; sys.argv = ['ltrtrace', '--version']\n"
            "from ltrtrace.cli import entry\n"
            "try:\n    entry()\nexcept SystemExit:\n    pass\n"
            "import sys; print(os.environ.get('NUMPY_MADVISE_HUGEPAGE'), 'numpy' in sys.modules)")
    for value, want in ((None, "0"), ("1", "1")):
        run_env = env if value is None else {**env, "NUMPY_MADVISE_HUGEPAGE": value}
        out = subprocess.run([sys.executable, "-c", code], env=run_env, capture_output=True, text=True, check=True)
        assert out.stdout.strip().splitlines()[-1] == f"{want} False"


def test_ci_checks_the_numbers():
    """The CI workflow lints, runs the fast tests and runs the golden pins (tests/test_golden.py, marked slow) with
    every marker, so a change of any pinned number fails CI."""
    from pathlib import Path

    workflow = Path(__file__).resolve().parents[1] / ".github" / "workflows" / "tests.yml"
    runs = [ln.split("run:", 1)[1].strip() for ln in workflow.read_text().splitlines() if ln.strip().startswith("run:")]
    assert "ruff check src tests" in runs and 'pytest -m "not slow"' in runs
    assert "pytest tests/test_golden.py" in runs, runs
