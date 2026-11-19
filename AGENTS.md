# Contributing to PydSEAMSlib

PydSEAMSlib builds the PyPI distribution `pydseamslib`: the Python
bindings of the d-SEAMS engine. This file is the contract for changes to
this repository.

## Layout

| Path | Contents |
| --- | --- |
| `src/pydseams/`, `src/bindings.cpp` | The `Frame` package and its nanobind module `pydseams.yoda`, which links `libyodaLib` |
| `src/pydseamslib/` | Re-exports `pydseams` |
| `src/dseams/` | `import dseams`: the nanobind module `_core` over `libseams` (`seams.h`), `frames.py` and `userframe.py` over `libseams_lang` |
| `src/seams/` | The ctypes package over `libseams`; not installed by the wheel |
| `examples/python/`, `examples/ctypes/` | Programs the docs include and `meson test` runs |
| `tests/python/` | Tests of `pydseams` and `pydseamslib` |
| `tests/dseams/` | Tests of `dseams`, the ctypes wrappers, the examples and the documentation pages |
| `docs/orgmode/` | Documentation source; `docs/export.el` writes the reStructuredText under `docs/source/` |
| `docs/source/` | Sphinx tree; the `.md` pages there are written by hand |
| `subprojects/` | Meson wraps: `seams-core` and the libraries it pins |
| `changelog.d/` | towncrier fragments |

## The engine

The C++ engine, its C API and the CLI live in
[d-SEAMS/seams-core](https://github.com/d-SEAMS/seams-core).
`subprojects/seams-core.wrap` pins one engine revision, the tag
`v2.14.1`.
The other wraps (`linkcell`, `minimage`, `readcon-*`, `vesin`, `gpulite`,
`eigen`) carry the revisions and patches seams-core pins at that revision.

Engine changes go to d-SEAMS/seams-core as their own pull request. Do
not patch `subprojects/seams-core/` here, and do not add a second engine
to this repository. Moving the pin is its own commit with a
`changed` fragment that names the new revision.

## Build and test

The engine links Rust crates, so `cargo` must be on `PATH`. Building the
wheel compiles the whole engine; use a build host, not a laptop.

```bash
# Development tree and the tests that read the build directory
meson setup bbdir --force-fallback-for=eigen3
meson compile -C bbdir
meson test -C bbdir --suite pydseamslib --print-errorlogs

# Wheel and sdist
python -m pip wheel --no-deps -w dist .
python -m build --sdist

# Installed package
python -m pip install dist/pydseamslib-*.whl pytest hypothesis ase
python -m pytest tests/python
python -m pytest --import-mode=importlib tests/dseams/test_dseams.py \
  tests/dseams/installed_package.py tests/dseams/docs tests/dseams/test_doc_snippets.py
python scripts/check_installed_api.py
```

`tests/python` imports sibling test modules, so it runs in pytest's
default import mode; `tests/dseams` runs with `--import-mode=importlib`.
Tests that read seams-core inputs take them from `SEAMS_CORE_SRC`, else
`subprojects/seams-core`.

`tests/dseams/frame_parity.py` compares `dseams.userframe` with the
published pydseamslib 2.11.1. Install 2.11.1 into a separate environment
and set `PYDSEAMS211` to its `site-packages` and `PYDSEAMS211_PYTHON` to
its interpreter.

A change that alters a computed number changes the code, or moves the
engine pin with a fragment that names the new value. Do not skip,
weaken, or delete a test to make a run pass.

## Formatting

`pre-commit` runs clang-format, `meson format`, and ruff. Commit a
formatter-only diff on its own with a `style:` subject.

## Commits and changelog

- Conventional Commits: `feat:`, `fix:`, `docs:`, `test:`, `build:`,
  `ci:`, `refactor:`, `style:`, `chore:`, with an optional scope.
  One logical change per commit.
- Every user-visible change adds a fragment in the same commit:
  `changelog.d/+slug.<type>.md` or `<issue>.<type>.md`, with `type` one
  of `added`, `changed`, `deprecated`, `removed`, `fixed`, `security`.
- The release cut is `towncrier build --version X.Y.Z`.

## Versions and release

The version is the seams-core release the wheel builds. It appears in
`meson.build`, `pyproject.toml`, `src/pydseams/__init__.py`,
`tests/python/test_bindings.py` and `CITATION.cff`; change them together.

1. d-SEAMS/seams-core tags `vX.Y.Z`.
2. `subprojects/seams-core.wrap` moves to that tag.
3. Cut the changelog, tag `vX.Y.Z` here.
4. `.github/workflows/build_wheels.yml` builds the sdist and the wheels
   and publishes them to PyPI through trusted publishing from
   d-SEAMS/PydSEAMSlib. The tag must equal `v` plus the
   `pyproject.toml` version.
