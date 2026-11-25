# Copyright (c) d-SEAMS developers
"""The compiled module's public names match the bindings and the docs."""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

from pydseams import yoda

ROOT = Path(__file__).resolve().parents[2]
PAGE = ROOT / "docs" / "source" / "api.md"


def _public_names(*, features=None):
    path = ROOT / "scripts" / "yoda_public_names.py"
    spec = importlib.util.spec_from_file_location("yoda_public_names", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.public_names(ROOT, features=features)


def _table_names(text: str) -> set[str]:
    section = text.split("## Compiled module", 1)[1]
    names: set[str] = set()
    for line in section.splitlines():
        if not line.startswith("| `"):
            continue
        cell = line.split("|", 2)[1]
        names.update(re.findall(r"`([^`]+)`", cell))
    return names


def test_package_exposes_the_bindings_this_build_registers():
    features = yoda._build_features
    names = _public_names(features=features)
    live = {name for name in dir(yoda) if not name.startswith("_")}
    n_optional = (
        int(features["with_chemfiles"])
        + int(features["with_readcon"])
        + 2 * int(features["with_ira"])
    )
    assert live == set(names.exposed)
    assert len(live) == 152 + n_optional
    assert len(names.exposed_functions) == 124 + n_optional
    assert len(names.exposed_types) == 26
    assert names.exposed_aliases == ("Kind", "Family")


def test_compiled_module_page_lists_every_registered_name():
    names = _public_names()
    text = PAGE.read_text(encoding="utf-8")
    assert _table_names(text) == set(names.registered)
    assert len(names.exposed) == 154
    assert "{{ yoda_names }}" in text
    assert "{{ yoda_functions }}" in text
    assert "{{ yoda_types }}" in text
    assert "134 names" not in text
