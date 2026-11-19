"""Release version for cogapp snippets.

``[project].version`` in ``pyproject.toml`` is the only place a person
writes the version. Every copy listed in ``ci/version/files.txt`` is
regenerated from it. Call cogapp as ``python -m cogapp``: cocogitto
also installs a ``cog`` executable.

    python ci/version/pydseams_version.py sync X.Y.Z
    python -m cogapp -r -I ci/version @ci/version/files.txt
"""

from __future__ import annotations

import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

with (ROOT / "pyproject.toml").open("rb") as fh:
    VERSION = tomllib.load(fh)["project"]["version"]


def sync(version: str) -> None:
    """Write ``version`` into ``pyproject.toml`` and nowhere else."""
    path = ROOT / "pyproject.toml"
    text = path.read_text(encoding="utf-8")
    new, count = re.subn(
        r'(?m)^version = "[^"]+"\n',
        f'version = "{version}"\n',
        text,
        count=1,
    )
    if count != 1:
        raise SystemExit("pyproject.toml has no project version key")
    path.write_text(new, encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "sync":
        sync(sys.argv[2])
    else:
        raise SystemExit(f"usage: {sys.argv[0]} sync VERSION")
