# Copyright (c) d-SEAMS developers
"""Public names registered by ``src/bindings.cpp``.

The tables on the compiled-module page and the documentation build read
these sets. The documentation counts refer to the published wheels.
"""

from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass
from pathlib import Path

# Binding name -> seams-core option that has to be enabled.
_GATED = {
    "readChemfiles": "with_chemfiles",
    "readCon": "with_readcon",
    "ira_match": "with_ira",
    "sofi_point_group": "with_ira",
}

_FUNCTION = re.compile(r'\bm\.def\(\s*"([^"]+)"')
_TYPE = re.compile(
    r'\bnb::(?:class_|enum_)<[^;]*?\(\s*m\s*,\s*"([^"]+)"',
    re.DOTALL,
)
_ALIAS = re.compile(r'\bm\.attr\(\s*"([^"]+)"\s*\)\s*=')
_SUBPROJECT = re.compile(
    r"subproject\(\s*'seams-core'.*?default_options:\s*\[(.*?)\]",
    re.DOTALL,
)
_OPTION = re.compile(r"^(?:-D|')(with_[A-Za-z0-9_]+)=(enabled|disabled|true|false)'?$")


@dataclass(frozen=True)
class YodaNames:
    """Unique public names, split the way the compiled-module page describes them."""

    functions: tuple[str, ...]
    types: tuple[str, ...]
    aliases: tuple[str, ...]
    exposed_functions: tuple[str, ...]
    exposed_types: tuple[str, ...]
    exposed_aliases: tuple[str, ...]

    @property
    def registered(self) -> tuple[str, ...]:
        return self.functions + self.types + self.aliases

    @property
    def exposed(self) -> tuple[str, ...]:
        return self.exposed_functions + self.exposed_types + self.exposed_aliases


def _unique(names: list[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    ordered: list[str] = []
    for name in names:
        if not name.startswith("_") and name not in seen:
            seen.add(name)
            ordered.append(name)
    return tuple(ordered)


def wheel_features(root: Path) -> dict[str, bool]:
    """Resolve literal features from the wheel and engine build settings."""
    meson_text = (root / "meson.build").read_text(encoding="utf-8")
    found = _SUBPROJECT.search(meson_text)
    if found is None:
        raise ValueError("meson.build has no seams-core subproject() options")
    wheel = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    setup_args = wheel["tool"]["cibuildwheel"]["config-settings"]["setup-args"]
    settings = [entry.strip() for entry in found.group(1).split(",")]
    options = dict(
        match.groups()
        for setting in [*settings, *setup_args]
        if (match := _OPTION.fullmatch(setting))
    )
    missing = sorted(set(_GATED.values()) - set(options))
    if missing:
        raise ValueError(f"seams-core subproject options missing {missing}")
    return {option: value in {"enabled", "true"} for option, value in options.items()}


def _exposed(names: tuple[str, ...], features: dict[str, bool]) -> tuple[str, ...]:
    return tuple(name for name in names if name not in _GATED or features[_GATED[name]])


def public_names(
    root: Path | None = None, *, features: dict[str, bool] | None = None
) -> YodaNames:
    """Binding names and the subset compiled with the specified features."""
    root = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    bindings = (root / "src" / "bindings.cpp").read_text(encoding="utf-8")
    if features is None:
        features = wheel_features(root)
    functions = _unique(_FUNCTION.findall(bindings))
    types = _unique(_TYPE.findall(bindings))
    aliases = _unique(_ALIAS.findall(bindings))
    return YodaNames(
        functions=functions,
        types=types,
        aliases=aliases,
        exposed_functions=_exposed(functions, features),
        exposed_types=_exposed(types, features),
        exposed_aliases=_exposed(aliases, features),
    )
