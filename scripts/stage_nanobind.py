# Copyright (c) d-SEAMS developers
"""Copy installed nanobind sources into the build directory."""

import pathlib
import shutil
import sys

src = pathlib.Path(sys.argv[1])
out = pathlib.Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
for name in sys.argv[3:]:
    shutil.copyfile(src / name, out / name)
