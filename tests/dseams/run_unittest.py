# Copyright (c) d-SEAMS developers
"""Import the built nanobind module and run unittest discovery."""

import argparse
import ast
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


def publish(source, dest):
    """Put source at dest without ever leaving a partial file there.

    meson test runs python_dseams and python_dseams_no_checkout at the same
    time, and both publish into this directory. Copying over dest rewrites
    the inode another process may have mapped or be loading, and that process
    dies with signal 11 or 7. A rename leaves the old inode intact for the
    processes that have it open.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(f".{dest.name}.{os.getpid()}.partial")
    try:
        shutil.copy(source, partial)
        os.replace(partial, dest)
    finally:
        partial.unlink(missing_ok=True)


pkg = stage / "dseams"
package = src / "src" / "dseams"
for name in ("__init__.py", "__main__.py", "frames.py", "userframe.py"):
    publish(package / name, pkg / name)
for so in cores:
    publish(so, pkg / so.name)

# Doc snippets import pydseams.yoda. The extension is built at the top of
# the build tree; put it next to a copy of the package so that import works.
front = stage / "pydseams"
source_front = src / "src" / "pydseams"
for path in source_front.rglob("*.py"):
    publish(path, front / path.relative_to(source_front))
for pattern in ("yoda*.so", "yoda*.pyd", "yoda*.dylib"):
    for shared in build.glob(pattern):
        publish(shared, front / shared.name)

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
    # The staged yoda extension is a copy, so its build rpath no longer
    # points at the vendored chemfiles library.
    chemfiles = engine / "chemfiles-prefix" / "lib"
    if chemfiles.is_dir():
        libdirs.append(str(chemfiles))
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


def _testcase_aliases(tree):
    names = set()
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module == "unittest":
            for alias in node.names:
                if alias.name == "TestCase":
                    names.add(alias.asname or alias.name)
    return names


def _is_unittest_testcase(class_node, aliases):
    for base in class_node.bases:
        if (
            isinstance(base, ast.Attribute)
            and base.attr == "TestCase"
            and isinstance(base.value, ast.Name)
            and base.value.id == "unittest"
        ):
            return True
        if isinstance(base, ast.Name) and base.id in aliases:
            return True
    return False


def _has_main_guard(tree):
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        test = node.test
        if (
            isinstance(test, ast.Compare)
            and len(test.ops) == 1
            and isinstance(test.ops[0], ast.Eq)
            and len(test.comparators) == 1
            and isinstance(test.comparators[0], ast.Constant)
            and test.comparators[0].value == "__main__"
        ):
            left = test.left
            if isinstance(left, ast.Name) and left.id == "__name__":
                return True
    return False


def module_ast_features(source):
    tree = ast.parse(source)
    aliases = _testcase_aliases(tree)
    pytest_like = False
    for node in tree.body:
        if isinstance(
            node, (ast.FunctionDef, ast.AsyncFunctionDef)
        ) and node.name.startswith("test_"):
            pytest_like = True
            break
        if (
            isinstance(node, ast.ClassDef)
            and node.name.startswith("Test")
            and not _is_unittest_testcase(node, aliases)
        ):
            pytest_like = True
            break
    script_only = _has_main_guard(tree) and not pytest_like
    return pytest_like, script_only


def resolve_module_path(name):
    path = root / module_filename(name)
    if path.is_file():
        return path
    return None


def should_route_pytest(path):
    source = path.read_text(encoding="utf-8")
    if "import pytest" in source:
        return True
    try:
        pytest_like, _ = module_ast_features(source)
    except SyntaxError:
        return False
    return pytest_like


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
failed = False
pytest_paths = []
for name in sorted(chosen - unittest_modules - failed_modules):
    path = resolve_module_path(name)
    if not path or not path.is_file():
        print(f"selected module file not found: {name}", file=sys.stderr)
        failed = True
        continue
    if should_route_pytest(path):
        pytest_paths.append(path)

missing_include = sorted(include - modules)
missing_exclude = sorted(exclude - modules)
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
