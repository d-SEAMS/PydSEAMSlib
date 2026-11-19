# Copyright (c) d-SEAMS developers
"""The dseams reference page matches the live call signatures.

``docs/source/reference/dseams.md`` is the hand-written symbol reference
(it is not an Org export). Each ``###`` heading inside the symbol block
is resolved on the imported ``dseams`` module. Parameter names, kinds,
defaults, and annotations are compared with that object.

Python callables go through ``inspect.signature``. Nanobind methods do
not implement it; their ``__nb_signature__`` (published by the extension
built into this process) is turned into an ``inspect.Signature``. A
documentation build stub (``DSEAMS_DOC_BUILD``) is refused.
"""

from __future__ import annotations

import ast
import inspect
import os
import re
import sys
import unittest
from pathlib import Path

os.environ.pop("DSEAMS_DOC_BUILD", None)

import dseams
from dseams import userframe

ROOT = Path(__file__).resolve().parents[2]
PAGE = ROOT / "docs" / "source" / "reference" / "dseams.md"
BLOCK = re.compile(
    r"<!-- API_SYMBOLS_BEGIN -->(.*)<!-- API_SYMBOLS_END -->",
    re.DOTALL,
)
HEADING = re.compile(r"^###\s+(.+)$", re.MULTILINE)
PARAM = re.compile(r"\*\*Parameters\*\*:\s*(.+)")
RETURNS = re.compile(r"\*\*Returns\*\*:\s*(.+)")
TOKEN = re.compile(
    r"`(?P<token>[^`]+)`"
    r"(?:\s*\((?:`)?(?P<type>int|float|str|bool|tuple|None)(?:`)?\))?"
)
PARAM_NAME = re.compile(
    r"^(?P<stars>\*{0,2})(?P<name>[A-Za-z_][A-Za-z0-9_]*)(?:\s*=\s*(?P<default>.+))?$"
)
PURPOSE_CODE = re.compile(r"\(`(-?\d+)`\)")
ANN = {
    "int": int,
    "float": float,
    "str": str,
    "bool": bool,
    "tuple": tuple,
    "None": type(None),
}


def _annotation(text):
    text = text.strip()
    if not text or text == "inspect._empty":
        return inspect.Signature.empty
    if text in ANN:
        return ANN[text]
    return text


def _split_args(inside):
    parts = []
    buf = []
    depth = 0
    for char in inside:
        if char == "," and depth == 0:
            parts.append("".join(buf).strip())
            buf = []
            continue
        if char in "([{":
            depth += 1
        elif char in ")]}":
            depth -= 1
        buf.append(char)
    tail = "".join(buf).strip()
    if tail:
        parts.append(tail)
    return parts


def _signature_from_nanobind(obj):
    """``inspect.Signature`` for one nanobind callable.

    ``inspect.signature`` raises ``ValueError`` on ``nb_method``. The
    extension stores the same text and default tuple on
    ``__nb_signature__``.
    """
    text, _doc, defaults = obj.__nb_signature__[0]
    if "->" in text:
        head, ret = text.split("->", 1)
        return_annotation = _annotation(ret.strip())
    else:
        head = text
        return_annotation = inspect.Signature.empty
    inside = head[head.index("(") + 1 : head.rindex(")")]
    params = []
    for part in _split_args(inside):
        if not part:
            continue
        name, _, rest = part.partition(":")
        name = name.strip()
        annotation = inspect.Parameter.empty
        default = inspect.Parameter.empty
        rest = rest.strip()
        if rest:
            if "=" in rest:
                ann_text, raw_default = rest.split("=", 1)
                annotation = _annotation(ann_text.strip())
                raw_default = raw_default.strip()
                if raw_default.startswith("\\") and raw_default[1:].isdigit():
                    default = defaults[int(raw_default[1:])]
                else:
                    default = ast.literal_eval(raw_default)
            else:
                annotation = _annotation(rest)
        if name.startswith("**"):
            kind = inspect.Parameter.VAR_KEYWORD
            name = name[2:]
        elif name.startswith("*"):
            kind = inspect.Parameter.VAR_POSITIONAL
            name = name[1:]
        else:
            kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
        params.append(
            inspect.Parameter(name, kind, default=default, annotation=annotation)
        )
    return inspect.Signature(params, return_annotation=return_annotation)


def live_signature(obj):
    """Signature of ``obj``, or ``None`` when it is not a callable.

    Constants and dict subclasses that do not define ``__init__`` return
    ``None``. Properties return an empty signature.
    """
    if isinstance(obj, property):
        return inspect.Signature()
    if inspect.isclass(obj):
        init = obj.__dict__.get("__init__", obj.__init__)
        if hasattr(init, "__nb_signature__"):
            return _signature_from_nanobind(init)
        if "__init__" not in obj.__dict__:
            try:
                return inspect.signature(obj)
            except (ValueError, TypeError):
                return None
        return inspect.signature(obj)
    if hasattr(obj, "__nb_signature__"):
        try:
            inspect.signature(obj)
        except (ValueError, TypeError):
            return _signature_from_nanobind(obj)
        else:
            return inspect.signature(obj)
    if not callable(obj):
        return None
    return inspect.signature(obj)


def _lookup_name(name, obj):
    modules = []
    mod_name = getattr(obj, "__module__", None)
    if mod_name and mod_name in sys.modules:
        modules.append(sys.modules[mod_name])
    modules.extend([dseams, userframe])
    seen = set()
    for module in modules:
        ident = id(module)
        if ident in seen:
            continue
        seen.add(ident)
        if hasattr(module, name):
            return getattr(module, name)
    raise LookupError(name)


def _default_value(raw, obj):
    try:
        return ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return _lookup_name(raw, obj)


def documented_parameters(line, obj):
    """Parameter tuples from one ``**Parameters**`` line.

    Each item is ``(name, kind, default, annotation)``. ``self`` is kept
    so the comparison can drop it on both sides. A line whose text is
    ``none`` yields an empty list.
    """
    body = line.strip().rstrip(".")
    if body.lower() == "none":
        return []
    found = list(TOKEN.finditer(body))
    if not found:
        raise AssertionError(f"no parameters in {line!r}")
    params = []
    for match in found:
        token = match.group("token").strip()
        parsed = PARAM_NAME.match(token)
        if not parsed:
            # Shape notes such as `(n, 3)` sit in backticks beside the name.
            if re.search(r"[^A-Za-z0-9_*=]", token):
                continue
            raise AssertionError(f"unparsed parameter token {token!r} in {line!r}")
        type_name = match.group("type")
        annotation = ANN[type_name] if type_name else inspect.Parameter.empty
        name = parsed.group("name")
        stars = parsed.group("stars")
        raw_default = parsed.group("default")
        if stars == "**":
            kind = inspect.Parameter.VAR_KEYWORD
        elif stars == "*":
            kind = inspect.Parameter.VAR_POSITIONAL
        else:
            kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
        if raw_default is None:
            default = inspect.Parameter.empty
        else:
            default = _default_value(raw_default.strip(), obj)
        params.append((name, kind, default, annotation))
    return params


def _strip_self(params):
    if params and params[0][0] == "self":
        return params[1:]
    return params


def _live_parameters(sig):
    params = []
    for param in sig.parameters.values():
        params.append((param.name, param.kind, param.default, param.annotation))
    return _strip_self(params)


def resolve_heading(text):
    names = re.findall(r"`([^`]+)`", text)
    if not names:
        raise AssertionError(f"heading has no symbol: {text}")
    objects = []
    for name in names:
        parts = name.split(".")
        if parts[0] == "dseams":
            obj = dseams
            parts = parts[1:]
        elif parts[0] == "userframe":
            obj = userframe
            parts = parts[1:]
        elif hasattr(dseams, parts[0]):
            obj = dseams
        else:
            raise AssertionError(f"cannot resolve {name}")
        for part in parts:
            obj = getattr(obj, part)
        objects.append((name, obj))
    return objects


def _sections(text):
    match = BLOCK.search(text)
    if not match:
        raise AssertionError("API_SYMBOLS markers are missing")
    body = match.group(1)
    marks = list(HEADING.finditer(body))
    sections = []
    for index, mark in enumerate(marks):
        end = marks[index + 1].start() if index + 1 < len(marks) else len(body)
        sections.append((mark.group(1).strip(), body[mark.end() : end]))
    return sections


def _same_default(documented, live):
    if documented is inspect.Parameter.empty or live is inspect.Parameter.empty:
        return documented is live
    return documented == live and type(documented) is type(live)


class ReferenceSignatureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        core = dseams._core
        if not getattr(core, "__file__", None):
            raise AssertionError(
                "dseams._core is the documentation stub, not the extension"
            )
        # The extension must be the one this interpreter loaded.
        if not str(core.__file__).endswith((".so", ".pyd", ".dylib")):
            raise AssertionError(
                f"dseams._core is not a compiled module: {core.__file__}"
            )

    def test_every_documented_signature(self):
        text = PAGE.read_text(encoding="utf-8")
        failures = []
        seen = []
        for heading, body in _sections(text):
            param_match = PARAM.search(body)
            returns_match = RETURNS.search(body)
            for name, obj in resolve_heading(heading):
                seen.append(name)
                sig = live_signature(obj)
                if param_match is None:
                    if sig is not None and _live_parameters(sig):
                        failures.append(
                            f"{name}: documented no parameters, live {_live_parameters(sig)}"
                        )
                    continue
                try:
                    documented = _strip_self(
                        documented_parameters(param_match.group(1), obj)
                    )
                except (AssertionError, LookupError, SyntaxError, ValueError) as exc:
                    failures.append(f"{name}: {exc}")
                    continue
                if sig is None:
                    if documented:
                        failures.append(
                            f"{name}: not callable but documents {documented}"
                        )
                    if isinstance(obj, int):
                        codes = PURPOSE_CODE.findall(body)
                        if not codes or int(codes[0]) != obj:
                            failures.append(
                                f"{name}: purpose integer {codes} != live {obj}"
                            )
                    continue
                live = _live_parameters(sig)
                if len(documented) != len(live):
                    failures.append(
                        f"{name}: docs {[item[0] for item in documented]} "
                        f"!= live {[item[0] for item in live]}"
                    )
                for doc_param, live_param in zip(documented, live):
                    d_name, d_kind, d_default, d_ann = doc_param
                    l_name, l_kind, l_default, l_ann = live_param
                    if (
                        d_name != l_name
                        or d_kind != l_kind
                        or not _same_default(d_default, l_default)
                    ):
                        failures.append(
                            f"{name}: docs ({d_name}, {d_kind}, {d_default!r}) "
                            f"!= live ({l_name}, {l_kind}, {l_default!r})"
                        )
                    if l_ann is not inspect.Parameter.empty:
                        if d_ann is inspect.Parameter.empty:
                            failures.append(
                                f"{name}: {l_name} is {l_ann.__name__} in the live signature"
                            )
                        elif d_ann is not l_ann:
                            failures.append(
                                f"{name}: {l_name} documented {d_ann} != live {l_ann}"
                            )
                ret = sig.return_annotation
                if (
                    ret is not inspect.Signature.empty
                    and ret is not None
                    and ret is not type(None)
                    and returns_match is not None
                ):
                    token = f"`{ret.__name__}`"
                    if token not in returns_match.group(1):
                        failures.append(
                            f"{name}: return annotation {ret.__name__} missing from {returns_match.group(1)!r}"
                        )
        self.assertGreater(len(seen), 30)
        # Python functions must have been read with inspect.signature.
        # A nanobind-only path would miss dseams.cages, which is pure Python.
        cages_sig = inspect.signature(dseams.cages)
        self.assertEqual(
            [param.name for param in cages_sig.parameters.values()],
            ["positions", "box", "cutoff", "graph", "k"],
        )
        self.assertEqual(cages_sig.parameters["graph"].default, "seeded")
        self.assertNotIn("ring_adjacent", cages_sig.parameters)
        if failures:
            self.fail("\n".join(failures))


if __name__ == "__main__":
    unittest.main()
