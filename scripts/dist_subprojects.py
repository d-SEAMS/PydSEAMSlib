# Copyright (c) d-SEAMS developers
"""Carry patched wrap-git subprojects into a `meson dist` tree.

`meson dist --include-subprojects` archives the committed tree of a
wrap-git subproject, which drops the wrap's diff_files patches and any
packagefile overlay. linkcell needs its patches (lc_mesh_cut,
lc_partition_weights) for seams-core to compile. This copies each
subproject's modified and untracked files from the configured source
tree, removes the files a patch deletes, and drops VCS metadata.
"""

import os
import shutil
import subprocess
from pathlib import Path

source = Path(os.environ["MESON_SOURCE_ROOT"]) / "subprojects"
dist = Path(os.environ["MESON_DIST_ROOT"]) / "subprojects"


def _git(repo, *args):
    out = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return [line for line in out.split("\0") if line]


if dist.is_dir():
    for child in sorted(dist.iterdir()):
        src = source / child.name
        if child.is_dir() and (src / ".git").exists():
            changed = _git(src, "ls-files", "-z", "--modified", "--others", "--exclude-standard")
            deleted = set(_git(src, "ls-files", "-z", "--deleted"))
            for rel in changed:
                if rel in deleted:
                    continue
                dest = child / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src / rel, dest)
            for rel in deleted:
                (child / rel).unlink(missing_ok=True)
        vcs = child / ".git"
        if vcs.is_dir():
            shutil.rmtree(vcs)
        elif vcs.exists():
            vcs.unlink()
