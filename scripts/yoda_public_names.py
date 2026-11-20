# Copyright (c) d-SEAMS developers
"""Public names registered by ``src/bindings.cpp``.

The tables on the compiled-module page and the documentation build read
these sets. A name inside a compile gate is part of this package's build
only when ``meson.build`` enables that seams-core option.
"""

from __future__ import annotations

import re
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
_OPTION = re.compile(r"'(with_[A-Za-z0-9_]+)=(enabled|disabled|true|false)'")


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
        if name not in seen:
            seen.add(name)
            ordered.append(name)
    return tuple(ordered)


def subproject_options(meson_text: str) -> dict[str, str]:
    """``with_*`` values from the seams-core ``subproject()`` call."""
    found = _SUBPROJECT.search(meson_text)
    if found is None:
        raise ValueError("meson.build has no seams-core subproject() options")
    options = dict(_OPTION.findall(found.group(1)))
    missing = sorted(set(_GATED.values()) - set(options))
    if missing:
        raise ValueError(f"seams-core subproject options missing {missing}")
    return options


def _enabled(option: str, options: dict[str, str]) -> bool:
    return options[option] in {"enabled", "true"}


def _exposed(names: tuple[str, ...], options: dict[str, str]) -> tuple[str, ...]:
    return tuple(
        name
        for name in names
        if _GATED.get(name) is None or _enabled(_GATED[name], options)
    )


def public_names(root: Path | None = None) -> YodaNames:
    """Names in ``src/bindings.cpp``, and the subset this package builds."""
    root = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    bindings = (root / "src" / "bindings.cpp").read_text(encoding="utf-8")
    options = subproject_options((root / "meson.build").read_text(encoding="utf-8"))
    functions = _unique(_FUNCTION.findall(bindings))
    types = _unique(_TYPE.findall(bindings))
    aliases = _unique(_ALIAS.findall(bindings))
    return YodaNames(
        functions=functions,
        types=types,
        aliases=aliases,
        exposed_functions=_exposed(functions, options),
        exposed_types=_exposed(types, options),
        exposed_aliases=_exposed(aliases, options),
    )
