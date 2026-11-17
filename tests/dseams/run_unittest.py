# Copyright (c) d-SEAMS developers
"""Import the built nanobind module and run unittest discovery."""

import argparse
import fnmatch
import os
import re
import shutil
import sys
import unittest
from pathlib import Path

try:
    import numpy  # noqa: F401
except ImportError:
    print("numpy is not installed")
    sys.exit(77)

build = Path(os.environ.get("DSEAMS_BUILD", "build"))
src = Path(os.environ.get("DSEAMS_SRC", Path(__file__).resolve().parents[2]))
# A conda python and the extension must share a stable ABI tag. A module
# built for another interpreter is a skip, not a failed import.
tag = sys.implementation.cache_tag  # cpython-312
# With python.allow_limited_api the module is _core.abi3.so, which any
# CPython from 3.12 loads.
abi3 = sys.version_info >= (3, 12)
cores = [
    path
    for path in (build / "src" / "dseams").glob("_core*.so")
    if tag in path.name or (abi3 and ".abi3." in path.name)
]
if not cores:
    built = [path.name for path in build.rglob("_core*.so")]
    print(f"dseams._core was not built for {tag}; found {built or 'nothing'}")
    sys.exit(77)

stage = build / "python-stage"
pkg = stage / "dseams"
pkg.mkdir(parents=True, exist_ok=True)
package = src / "src" / "dseams"
for name in ("__init__.py", "__main__.py", "frames.py", "userframe.py"):
    shutil.copy(package / name, pkg / name)
for so in cores:
    shutil.copy(so, pkg / so.name)

# Doc snippets import pydseams.yoda. The extension is built at the top of
# the build tree; put it next to a copy of the package so that import works.
front = stage / "pydseams"
front.mkdir(parents=True, exist_ok=True)
source_front = src / "src" / "pydseams"
for path in source_front.rglob("*.py"):
    dest = front / path.relative_to(source_front)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(path, dest)
for pattern in ("yoda*.so", "yoda*.pyd", "yoda*.dylib"):
    for shared in build.glob(pattern):
        shutil.copy(shared, front / shared.name)

if os.environ.get("SEAMS_LIB") and not os.environ.get("SEAMS_LIBDIR"):
    os.environ["SEAMS_LIBDIR"] = str(Path(os.environ["SEAMS_LIB"]).resolve().parent)

if os.environ.get("DSEAMS_PY_READY") != "1":
    env = os.environ.copy()
    libdirs = []
    for name in ("SEAMS_LIB", "SEAMS_LANG_LIB", "SEAMS_CORE_LIB"):
        value = env.get(name)
        if value:
            libdirs.append(str(Path(value).resolve().parent))
    engine = build / "subprojects" / "seams-core"
    libdirs.extend([str(engine / "src"), str(engine / "bindings")])
    libs = []
    for libdir in libdirs:
        if libdir and libdir not in libs:
            libs.append(libdir)
    env["LD_LIBRARY_PATH"] = os.pathsep.join(libs + [env.get("LD_LIBRARY_PATH", "")])
    env["PYTHONPATH"] = os.pathsep.join(
        [str(stage), str(src / "src"), env.get("PYTHONPATH", "")]
    )
    env["DSEAMS_PY_READY"] = "1"
    os.execve(sys.executable, [sys.executable, __file__, *sys.argv[1:]], env)

parser = argparse.ArgumentParser()
parser.add_argument("--include", nargs="*", default=[])
parser.add_argument("--exclude", nargs="*", default=[])
args = parser.parse_args()


def module_filename(name):
    return name if name.endswith(".py") else f"{name}.py"


include = {module_filename(name) for name in args.include}
exclude = {module_filename(name) for name in args.exclude}

# Same files unittest discovery visits: test_*.py, and packages below them.
_module_file = re.compile(r"[_a-z]\w*\.py$", re.IGNORECASE)


def discovery_modules(start):
    found = set()
    for path in sorted(start.iterdir()):
        if path.is_file():
            if _module_file.match(path.name) and fnmatch.fnmatch(
                path.name, "test_*.py"
            ):
                found.add(path.name)
        elif (path / "__init__.py").is_file():
            found.update(discovery_modules(path))
    return found


root = src / "tests" / "dseams"
modules = discovery_modules(root)
suite = unittest.defaultTestLoader.discover(str(root), pattern="test_*.py")


def iter_cases(node):
    for item in node:
        if isinstance(item, unittest.TestSuite):
            yield from iter_cases(item)
        else:
            yield item


def case_module(case):
    if isinstance(case, unittest.loader._FailedTest):
        stem = (case._testMethodName or "").split(".")[-1]
    else:
        stem = case.__class__.__module__.split(".")[-1]
    if not stem:
        return ""
    return module_filename(stem)


selected = unittest.TestSuite()
unittest_modules = set()
failed_modules = set()
for case in iter_cases(suite):
    name = case_module(case)
    if name:
        modules.add(name)
    # An import failure is not a passing module. Keep it no matter which
    # file the filters name, so a broken module cannot be selected away.
    if isinstance(case, unittest.loader._FailedTest):
        if name:
            failed_modules.add(name)
        selected.addTest(case)
        continue
    if name:
        unittest_modules.add(name)
    if include and name not in include:
        continue
    if name in exclude:
        continue
    selected.addTest(case)

chosen = set(modules)
if include:
    chosen &= include
chosen -= exclude
pytest_paths = []
for name in sorted(chosen - unittest_modules - failed_modules):
    path = root / name
    if path.is_file() and "import pytest" in path.read_text(encoding="utf-8"):
        pytest_paths.append(path)

missing_include = sorted(include - modules)
missing_exclude = sorted(exclude - modules)
failed = False
if missing_include:
    print(
        "include matched no test module: " + ", ".join(missing_include),
        file=sys.stderr,
    )
    failed = True
if missing_exclude:
    print(
        "exclude matched no test module: " + ", ".join(missing_exclude),
        file=sys.stderr,
    )
    failed = True
if (include or exclude) and selected.countTestCases() == 0 and not pytest_paths:
    print("filters selected no tests", file=sys.stderr)
    failed = True
else:
    if selected.countTestCases() or not (include or exclude):
        result = unittest.TextTestRunner(verbosity=1).run(selected)
        failed = failed or not result.wasSuccessful()
    if pytest_paths:
        import pytest

        status = int(
            pytest.main(
                [
                    "-q",
                    "-p",
                    "no:cacheprovider",
                    "--noconftest",
                    *[str(path) for path in pytest_paths],
                ]
            )
        )
        if status not in (0, 5):
            failed = True
        elif status == 5 and selected.countTestCases() == 0 and (include or exclude):
            print("filters selected no tests", file=sys.stderr)
            failed = True
sys.exit(1 if failed else 0)
