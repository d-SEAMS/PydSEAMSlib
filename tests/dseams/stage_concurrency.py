# Copyright (c) d-SEAMS developers
"""Concurrent run_unittest.py processes must not corrupt each other's stage.

run_unittest.py copies dseams._core and the Python sources into
<build>/python-stage and imports them from there. meson test starts
python_dseams and python_dseams_no_checkout at the same time, so two
processes publish the same files into the same directory. A copy that
rewrites the shared object in place corrupts the mapping of a process that
has already loaded it, or hands dlopen a half-written file, and the
interpreter dies with signal 11 or 7. The crash comes and goes between runs
because it depends on how the processes interleave, and a MALLOC_PERTURB_
value does not change that.

This starts many run_unittest.py processes at once, each with its own
MALLOC_PERTURB_, and requires every one to exit 0.
"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "run_unittest.py"
RUNNER_ARGV = [
    sys.executable,
    str(RUNNER),
    "--exclude",
    "test_seams_package.py",
    "test_doc_snippets.py",
]
PERTURB = ("0", "20", "139", "255")
PROCESSES = 12
ROUNDS = 10


def perturb(round_no, index):
    return PERTURB[(round_no + index) % len(PERTURB)]


def describe(returncode):
    if returncode < 0:
        return f"killed by signal {-returncode}"
    return f"exit status {returncode}"


def run_runner(env):
    return subprocess.run(
        RUNNER_ARGV,
        check=False,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=300,
    )


class StageConcurrency(unittest.TestCase):
    def test_published_files_are_replaced_not_rewritten(self):
        env = os.environ.copy()
        env.pop("DSEAMS_PY_READY", None)
        env["MALLOC_PERTURB_"] = "20"
        first = run_runner(env)
        self.assertEqual(0, first.returncode, first.stdout)
        build = Path(env["DSEAMS_BUILD"])
        shared = sorted((build / "python-stage" / "dseams").glob("_core*.so"))
        self.assertTrue(shared, "no dseams._core in the stage")
        before = {path: path.stat().st_ino for path in shared}
        # A process that has the shared object mapped keeps reading the
        # inode it opened. Rewriting that inode in place corrupts it.
        held = {path: path.open("rb") for path in shared}
        try:
            second = run_runner(env)
            self.assertEqual(0, second.returncode, second.stdout)
            for path, handle in held.items():
                self.assertNotEqual(
                    before[path],
                    path.stat().st_ino,
                    f"{path.name} was rewritten in place",
                )
                handle.seek(0, os.SEEK_END)
                self.assertGreater(handle.tell(), 0, f"{path.name} was truncated")
        finally:
            for handle in held.values():
                handle.close()

    def test_concurrent_runs_share_the_stage(self):
        env = os.environ.copy()
        env.pop("DSEAMS_PY_READY", None)
        failures = []
        for round_no in range(ROUNDS):
            children = []
            for index in range(PROCESSES):
                child_env = env.copy()
                child_env["MALLOC_PERTURB_"] = perturb(round_no, index)
                children.append(
                    subprocess.Popen(
                        RUNNER_ARGV,
                        env=child_env,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                    )
                )
            for index, child in enumerate(children):
                output, _ = child.communicate(timeout=300)
                if child.returncode != 0:
                    failures.append(
                        f"round {round_no} process {index} "
                        f"MALLOC_PERTURB_={perturb(round_no, index)}: "
                        f"{describe(child.returncode)}\n{output[-2000:]}"
                    )
        self.assertEqual([], failures, "\n".join(failures))


if __name__ == "__main__":
    unittest.main()
