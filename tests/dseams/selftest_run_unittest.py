# Copyright (c) d-SEAMS developers
"""Exercise run_unittest.py routing against a temporary stub tree."""

import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

if importlib.util.find_spec("pytest") is None:
    print("pytest is not installed", file=sys.stderr)
    sys.exit(77)

runner = Path(__file__).resolve().parent / "run_unittest.py"
repo = Path(__file__).resolve().parents[2]
package = repo / "src" / "dseams"
pydseams = repo / "src" / "pydseams"


def invoke(args, src_root, expect):
    env = os.environ.copy()
    env["DSEAMS_SRC"] = str(src_root)
    if "DSEAMS_BUILD" not in env:
        print("DSEAMS_BUILD is not set", file=sys.stderr)
        sys.exit(1)
    cmd = [sys.executable, str(runner), *args]
    proc = subprocess.run(cmd, check=False, env=env, capture_output=True, text=True)
    outline = proc.stdout + proc.stderr
    if proc.returncode != expect:
        print("unexpected exit", proc.returncode, "expected", expect, file=sys.stderr)
        print("args:", args, file=sys.stderr)
        print(outline, file=sys.stderr)
        sys.exit(1)


def _populate_src_root(src):
    stage_pkg = src / "src" / "dseams"
    stage_pkg.mkdir(parents=True)
    for name in ("__init__.py", "__main__.py", "frames.py", "userframe.py"):
        shutil.copy(package / name, stage_pkg / name)
    stage_front = src / "src" / "pydseams"
    for path in pydseams.rglob("*.py"):
        dest = stage_front / path.relative_to(pydseams)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(path, dest)


@contextmanager
def stub_tree(modules):
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp)
        _populate_src_root(src)
        tests = src / "tests" / "dseams"
        tests.mkdir(parents=True)
        for name, source in modules.items():
            path = tests / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(source, encoding="utf-8")
        yield src


def main():
    with stub_tree(
        {
            "test_plain_fail.py": "def test_always_fails():\n    assert False\n",
        }
    ) as src:
        invoke(["--include", "test_plain_fail.py"], src, 1)

    with stub_tree(
        {
            "test_plain_fail.py": "def test_plain():\n    assert False\n",
            "test_script_only.py": (
                'def main():\n    return\n\nif __name__ == "__main__":\n    main()\n'
            ),
        }
    ) as src:
        invoke(["--exclude", "test_script_only.py"], src, 1)

    with stub_tree(
        {
            "test_plain_pass.py": "def test_always_passes():\n    assert True\n",
        }
    ) as src:
        invoke(["--include", "test_plain_pass.py"], src, 0)

    with stub_tree({}) as src:
        invoke(["--include", "test_no_such_module.py"], src, 1)

    with stub_tree(
        {
            "pkg/__init__.py": "",
            "pkg/test_nested.py": "def test_nested_passes():\n    assert True\n",
        }
    ) as src:
        invoke(["--include", "test_nested.py"], src, 1)

    with stub_tree(
        {
            "test_plain_pass.py": "def test_always_passes():\n    assert True\n",
        }
    ) as src:
        invoke(
            ["--include", "test_plain_pass.py", "--exclude", "test_plain_pass.py"],
            src,
            1,
        )

    with stub_tree(
        {
            "test_broken_import.py": (
                "import not_a_real_stub_module\n\n"
                "def test_never_reached():\n    assert True\n"
            ),
        }
    ) as src:
        invoke(["--include", "test_broken_import.py"], src, 1)

    with stub_tree(
        {
            "test_script_only.py": (
                'def main():\n    return\n\nif __name__ == "__main__":\n    main()\n'
            ),
        }
    ) as src:
        invoke(["--include", "test_script_only.py"], src, 1)


if __name__ == "__main__":
    main()
