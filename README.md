# PydSEAMSlib

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/d-SEAMS/PydSEAMSlib/blob/main/notebooks/classify_ice.ipynb)

<p align="center">
  <img src="docs/source/_static/logo/pydseams-icon.png" alt="pydseams" width="96">
</p>

[![built with nix](https://builtwithnix.org/badge.svg)](https://builtwithnix.org)

Python bindings for the [d-SEAMS](https://dseams.info) C++ engine
([`seams-core`](https://github.com/d-SEAMS/seams-core)).

This repository builds the PyPI distribution `pydseamslib`. One wheel
holds three import names over the same engine:

| Import | What it is |
| --- | --- |
| `pydseams` | The `Frame` package: `read`, `from_ase`, `features`, `band`, solvis, and the compiled `pydseams.yoda` surface over `libyodaLib` |
| `pydseamslib` | Re-exports `pydseams` (`pydseamslib.Frame is pydseams.Frame`) |
| `dseams` | The 2.14 API over the stable C API (`seams.h`): `chill_plus`, `q4q6`, `neighbours`, `rings`, `cages`, `symmetry`, `dseams.Frame`, and `dseams.userframe` |

The C++ engine and `seams` CLI live in
[`seams-core`](https://github.com/d-SEAMS/seams-core), built here as a
Meson subproject; `libseams` and `libseams_lang` ship inside the wheel.
`src/seams` is the ctypes package over `libseams` (not in the wheel).
Lua/Fennel is `dseams` in [`yodaStruct`](https://github.com/d-SEAMS/yodaStruct).
Neighbour search is [`linkcell`](https://github.com/d-SEAMS/linkcell).
Do not grow a second engine here.

Released versions live in [CHANGELOG.md](CHANGELOG.md) (Keep a Changelog).
Unreleased notes are files under `changelog.d/`. The cut is
`towncrier build --version X.Y.Z`.

```bash
python -m pip install 'pydseamslib==2.14.1'
python -m pip install 'pydseamslib[ase]==2.14.1'      # ASE Atoms
python -m pip install 'pydseamslib[solvis]==2.14.1'   # solvis / PyVista
```

Each command installs a wheel containing the bundled engine and all three
import names: `dseams`, `pydseams`, and `pydseamslib`. No separate
seams-core installation is needed.

Nix flake:

```bash
nix build
nix develop
```

Docs: `docs/orgmode/` (ox-rst) and `docs/source/` (Shibuya).
Site: <https://python.dseams.info/>.

Book: [Python](https://docs.dseams.info/stable/python.html) and the
[Python Frame](https://docs.dseams.info/reference/python-frame.html)
reference in the d-SEAMS book.
Video: [pydseams walkthrough](https://odysee.com/dseams-python:be29b99d4059bf4afe45ad0a9ca71e20bf11beaa)
on Odysee.

```python
import pydseams as ds

frame = ds.read(
    "water.lammpstrj"
)  # release wheels also read .xyz, .con, .pdb, .gro, and .dcd
print(frame.chill_plus())
print(frame.cages())
print(frame.density(bins=100, axis="z"))

mixed = ds.read("ions.lammpstrj", all_atoms=True, atom_type=1)
sites = ds.yoda.parseSiteSpec("1=polar,2=apolar")
print(mixed.domain(sites, ds.yoda.Kind.polar))

ions = ds.yoda.parseSiteSpec("1=cationHead,2=anion")
print(mixed.pairs(ions))

frame = ds.from_ase(atoms)
atoms = frame.to_ase()
system = frame.to_solvis()  # optional extra
```

`ds.read` picks the engine reader from the suffix. `yoda` is the compiled
surface. Helpers (`Frame`, `io`, ASE, solvis) stay in Python.
`_core` and `cyoda` still name the same module.
`ds.available_readers()` reports the readers compiled into a wheel. The
2.14.1 PyPI wheels link chemfiles and readcon-core, so `.pdb`, `.gro`,
`.dcd`, and eOn `.con` all read directly next to LAMMPS and `.xyz`. Only
a source build left on the default `auto` detection, with one of those
libraries unavailable, raises `RuntimeError` for the suffixes that
library would have read.

`Frame.density`, `Frame.pairs`, and `Frame.domain` expose the same
site-resolved density, ion-pair, and connected-domain analyses as the CLI.
Use `ds.read(..., all_atoms=True)` when a site analysis needs every LAMMPS
type. A file with more than one element symbol gives each symbol its own
`c_type`, and the cutoff graph contains like and unlike pairs. Pass
`periodic=False` for a molecule so pairs are not wrapped across the
bounding span. `atom_type` is the species single-type ice routines still
receive. `pydseams.band.band_rows` names the pairs that enter or leave
that graph from one frame to the next.
ASE adapters accept nonsingular cells periodic in all three directions and
preserve cell orientation and displacement on roundtrip. Hydrogen donors use
an ASE ``mol-id`` array when present and periodic nearest-atom ownership
otherwise. All-atom ASE imports (``select=None``) use cutoff bonding;
hydrogen-bond topology requires a heavy-atom selection such as ``select="O"``.

Cutoff, frame, and *k* follow the same twelve-factor table as
`seams`: `SEAMS_CONFIG` or `./seams.env`, then the environment, then
the function argument. `pydseams.config` is the reader.

```python
import dseams

a = 2.75 * 4 / 3**0.5
sites = (
    (0, 0, 0),
    (0.5, 0.5, 0),
    (0.5, 0, 0.5),
    (0, 0.5, 0.5),
    (0.25, 0.25, 0.25),
    (0.75, 0.75, 0.25),
    (0.75, 0.25, 0.75),
    (0.25, 0.75, 0.75),
)
xyz = [[x * a, y * a, z * a] for x, y, z in sites]
assert dseams.chill_plus(xyz, [a, a, a])[1][0] == 8
```

`dseams` reads C-contiguous NumPy arrays in place. `box` is three edge
lengths or the nine LAMMPS numbers `ox, oy, oz, lx, ly, lz, xy, xz, yz`.
`cages`, `signature_cages`, and `symmetry` go through `libseams_lang`.

## Upgrading from 2.11

- `import pydseams` and `import pydseamslib` keep the 2.11 API, including
  `pydseams.yoda`, `cyoda`, and `_core`.
- CHILL+ on a cutoff list follows the seams-core 2.14 shell-count fix. On
  `tests/data/exampleTraj.lammpstrj` (oxygen, type 2) `Frame.chill_plus()`
  reports 9 interClathrate and 241 water, where 2.11 reported 12 and 238.
- `FrameFingerprint` and `KeyLibrary` also expose `colored`; `coloured`
  still works.
- `dseams` is new. `dseams.userframe.Frame` is a second `Frame` over
  `libseams_lang` with a `graph` argument on `chill_plus`; it is not
  `pydseams.Frame`.

## Development

The engine is [seams-core](https://github.com/d-SEAMS/seams-core), a
Meson subproject here. `subprojects/seams-core.wrap` pins tag `v2.14.1`,
and the other wraps carry the revisions that tag pins. Engine changes
go to d-SEAMS/seams-core, not to this repository.
The contributing guide at the repository root has the layout, the test
commands and the commit conventions.

The build needs a C++20 compiler and `cargo` with Rust 1.88 or newer
(linkcell, minimage and readcon are Rust crates; readcon-core and
readcon-db declare `rust-version = "1.88"`). Installing the sdist
builds the same tree, so it needs them too:

```bash
meson setup bbdir --force-fallback-for=eigen3
meson compile -C bbdir
meson test -C bbdir --suite pydseamslib
```

To build against a local seams-core checkout, put it where the wrap
would: Meson uses an existing `subprojects/seams-core` directory and
does not fetch the pinned revision.

```bash
ln -s /path/to/seams-core subprojects/seams-core
meson setup --wipe bbdir --force-fallback-for=eigen3
```

The tests that read seams-core inputs (`input/traj/...`) take them from
`SEAMS_CORE_SRC`, else `subprojects/seams-core`.

The release path is
`docs/orgmode/contributing/developer/release.org`.

`llms.txt` at the repository root, also served at the docs site root,
lists the documentation pages.

Primary author: Ruhila S (https://orcid.org/0000-0002-5443-9356). The project started as PSF GSoC 2023 (`pyseams`).

Requires Python 3.11+. CPython 3.11 has a version-specific wheel. Python 3.12 and newer use
the CPython 3.12 stable ABI (one `abi3` wheel per platform). Free-threaded CPython has no
limited ABI and is not a target.

# License

[MIT](LICENSE).

Install the `readcon` extra to stream labelled trajectory frames. Coordinate-only
formats such as XTC and TRR require a topology with matching atom order.

```python
from pydseams import iter_trajectory

for frame in iter_trajectory(
    "water.xtc",
    topology="water.gro",
    residue_name="SOL",
    atom_name="OW",
    hydrogen_names=("HW1", "HW2"),
):
    print(frame.frame, len(frame.rings))
```
