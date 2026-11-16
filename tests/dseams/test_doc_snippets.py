# Copyright (c) d-SEAMS developers
"""Run the Python blocks of the documentation pages against the installed package.

Each listed org page's ``#+begin_src python`` blocks run as one script, in page
order, with the seams-core sources (``SEAMS_CORE_SRC``, else the
subproject checkout) as the working directory, so ``input/traj/...``
paths resolve. ``>>>`` blocks in reStructuredText pages run as doctests
with ``examples/python`` on the path. A page whose imports need an
optional package that is not installed is skipped with that name.
"""

from __future__ import annotations

import importlib.util
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
ENGINE = Path(os.environ.get("SEAMS_CORE_SRC", REPO / "subprojects" / "seams-core"))
ORG_SRC = re.compile(
    r"^#\+begin_src python[^\n]*\n(.*?)^#\+end_src",
    re.MULTILINE | re.DOTALL | re.IGNORECASE,
)
RST_DOCTEST = re.compile(r"^\.\. doctest::\s*\n((?:\n|   .*\n)+)", re.MULTILINE)
OPTIONAL = ("ase", "MDAnalysis", "solvis")

# The pages written against the seams-core inputs. The other org pages
# show calls on placeholder files (water.lammpstrj) and are not run.
ORG_PAGES = [
    "docs/orgmode/howto/ase-integration.org",
    "docs/orgmode/howto/solvis-dseams.org",
    "docs/orgmode/tutorials/bulk-ice.org",
    "docs/orgmode/tutorials/confined-water.org",
    "docs/orgmode/tutorials/ice-and-clathrate.org",
    "docs/orgmode/tutorials/order-parameters.org",
    "docs/orgmode/tutorials/paper-figures.org",
    "docs/orgmode/tutorials/rings.org",
]
RST_PAGES = sorted(
    path.relative_to(REPO).as_posix()
    for path in (REPO / "docs" / "source" / "lang").glob("*.rst")
    if RST_DOCTEST.search(path.read_text(encoding="utf-8"))
)


def _env():
    env = os.environ.copy()
    env["DSEAMS_ROOT"] = str(ENGINE)
    env["SEAMS_CORE_SRC"] = str(ENGINE)
    if "SEAMS_LIBDIR" not in env:
        # core_counts loads libseams from the directory the wheel installs it in.
        spec = importlib.util.find_spec("dseams")
        if spec is not None and spec.origin:
            libs = Path(spec.origin).resolve().parents[1] / ".pydseamslib.mesonpy.libs"
            if libs.is_dir():
                env["SEAMS_LIBDIR"] = str(libs)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO / "examples" / "python"), env.get("PYTHONPATH", "")]
    )
    return env


def _run(cmd, cwd):
    done = subprocess.run(
        cmd, cwd=cwd, env=_env(), text=True, capture_output=True, check=False
    )
    if done.returncode != 0:
        raise AssertionError((done.stdout + done.stderr)[-4000:])
    return done.stdout


def _needs(body):
    for name in OPTIONAL:
        if (
            re.search(rf"^\s*(import|from) {name}\b", body, re.MULTILINE)
            and importlib.util.find_spec(name) is None
        ):
            pytest.skip(f"{name} is not installed")
    # Frame.to_solvis imports solvis itself (the [solvis] extra).
    if ".to_solvis(" in body and importlib.util.find_spec("solvis") is None:
        pytest.skip("solvis is not installed")


@pytest.mark.parametrize("page", ORG_PAGES)
def test_org_page(page, tmp_path):
    text = (REPO / page).read_text(encoding="utf-8")
    body = "\n".join(match.group(1) for match in ORG_SRC.finditer(text))
    _needs(body)
    script = tmp_path / "snippet.py"
    script.write_text(body, encoding="utf-8")
    with tempfile.TemporaryDirectory() as work:
        # Pages that write output (rdf_2d, monolayer_rings) write into a
        # scratch directory; input/ there links to the engine tree.
        os.symlink(ENGINE / "input", Path(work) / "input")
        _run([sys.executable, str(script)], cwd=work)


@pytest.mark.parametrize("page", RST_PAGES)
def test_rst_doctest(page, tmp_path):
    text = (REPO / page).read_text(encoding="utf-8")
    for index, match in enumerate(RST_DOCTEST.finditer(text)):
        lines = [line.removeprefix("   ") for line in match.group(1).splitlines()]
        doc = tmp_path / f"block{index}.txt"
        doc.write_text("\n".join(lines).strip() + "\n", encoding="utf-8")
        _run([sys.executable, "-m", "doctest", "-v", str(doc)], cwd=tmp_path)


def test_lang_tasks_python(monkeypatch):
    """The task-python page program on ice Ih: 256 six-rings, 128 hexagonal."""
    env = _env()
    if "SEAMS_LANG_LIB" not in os.environ and "SEAMS_LIBDIR" in env:
        name = "libseams_lang.dylib" if sys.platform == "darwin" else "libseams_lang.so"
        monkeypatch.setenv("SEAMS_LANG_LIB", str(Path(env["SEAMS_LIBDIR"]) / name))
    if "SEAMS_LANG_LIB" not in os.environ:
        pytest.skip("libseams_lang was not found; set SEAMS_LANG_LIB")
    out = _run(
        [sys.executable, str(REPO / "examples" / "python" / "lang_tasks_python.py")],
        cwd=ENGINE,
    )
    assert out.strip() == "pairs 512 rings 0 0 0 0 0 0 256 chill 0 128 0 0 0 0 0 0 0"
