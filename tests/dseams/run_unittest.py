# Copyright (c) d-SEAMS developers
"""Import the built nanobind module and run its unittest. Exit 77 if it was not built."""

import os
import shutil
import sys
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

if os.environ.get("DSEAMS_PY_READY") != "1":
    env = os.environ.copy()
    engine = build / "subprojects" / "seams-core"
    libs = [str(engine / "src"), str(engine / "bindings")]
    env["LD_LIBRARY_PATH"] = os.pathsep.join(libs + [env.get("LD_LIBRARY_PATH", "")])
    env["PYTHONPATH"] = os.pathsep.join(
        [str(stage), str(src / "src"), env.get("PYTHONPATH", "")]
    )
    env["DSEAMS_PY_READY"] = "1"
    os.execve(sys.executable, [sys.executable, __file__], env)

import unittest

suite = unittest.defaultTestLoader.discover(str(src / "tests" / "dseams"), pattern="test_*.py")
result = unittest.TextTestRunner(verbosity=1).run(suite)
sys.exit(0 if result.wasSuccessful() else 1)
