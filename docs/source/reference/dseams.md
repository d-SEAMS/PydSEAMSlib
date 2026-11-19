# dseams and the C API

Three Python modules call the same d-SEAMS engine:

- `dseams` (and `dseams.userframe`) is the nanobind and C-ABI module
  included in `pip install pydseamslib`. It reads C-contiguous NumPy
  arrays in place.
- `pydseams` (`pydseamslib`) is the high-level `Frame` package in
  `pip install pydseamslib`.
- `seams` is a pure-Python `ctypes` wrapper over `seams.h` for source
  checkouts (add `src` to `PYTHONPATH` and set `SEAMS_LIB` to
  `libseams.so` or `libseams.dylib`). Results are frozen dataclasses,
  and failing calls raise `seams.SeamsError` with `status` and
  `message`. `seams.from_pydseams(frame)` copies `frame.positions` and
  `frame.box` into a `System`.

## Quickstart

### `seams` (source checkout ctypes wrapper)

```python
import math

import numpy as np
from seams import PUBLISHED, System, set_num_threads

edge = 2.75 * 4.0 / math.sqrt(3.0)
fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
basis = ((0, 0, 0), (0.25, 0.25, 0.25))
positions = np.array(
    [[(s + b) * edge for s, b in zip(shift, site)] for shift in fcc for site in basis]
)

set_num_threads(1)
with System() as sys:
    sys.set_ortho(positions, (edge, edge, edge))
    ice = sys.chill_plus(rule=PUBLISHED)
    order = sys.q4q6()
    print(int(ice.counts[0]), float(order.q6.mean()))
    assert int(ice.counts[0]) == 8
    assert abs(float(order.q6.mean()) - 0.628539) < 1e-6
```

`examples/python/seams_ctypes_quick.py` runs that cell against a local `libseams` shared library. The program prints:

```
version 2.14.0
n 8
cubic 8
q6 0.628539
degree 4
six 0
fraction 0.0
cn 4.571429
```

Mean `q6` is `0.628539` and each atom has four neighbours. An 8-atom
unit cell is narrower than twice the cutoff, so cage census calls are
refused on it; a 2x2x2 supercell has a six-ring fraction above 0.5.

### `dseams`

These calls run directly against `pip install pydseamslib`, in a box
three cells wide so that topological cage and primitive-ring calls are
not clipped by the minimum-image cutoff:

```python
import math

import dseams
import numpy as np

edge = 2.75 * 4.0 / math.sqrt(3.0)
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
positions = np.array(
    [
        [(i + x) * edge, (j + y) * edge, (k + z) * edge]
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for x, y, z in sites
    ]
)
box = (3 * edge, 3 * edge, 3 * edge)
labels, counts = dseams.chill_plus(positions, box)
q4, q6 = dseams.q4q6(positions, box)
offsets, cols = dseams.neighbours(positions, box)
sizes, symmetry, hist, fraction = dseams.rings(positions, box)
per_atom, (n_hex, n_cubic, n_water, n_both) = dseams.cages(positions, box)
n_cages, n_atoms = dseams.signature_cages(positions, box, spec="512")
shell = dseams.symmetry(positions, box)
assert len(labels) == 216 and int(counts.sum()) == 216
assert len(q4) == 216 and len(q6) == 216
assert len(offsets) == 217
assert n_hex + n_cubic + n_water + n_both == 216
assert (n_cages, n_atoms) == (0, 0)
assert len(shell) == 216
```

`dseams.cages` returns a per-atom code and `(n_hex, n_cubic, n_water, n_both)`.
`dseams.signature_cages` counts closed polyhedra whose ring census
matches `spec` and returns `(n_cages, n_atoms)`. `dseams.symmetry`
returns the shell label: `0` ice Ih, `1` ice Ic, `2` stacking-disordered,
and `-1` otherwise.

`box` accepts three lengths or nine LAMMPS numbers (`ox, oy, oz, lx, ly, lz, xy, xz, yz`).
`dseams.PUBLISHED` (`"paper"`) averages $q_{lm}$ over the four closest
neighbours; `dseams.OVITO` (`"ovito"`) averages $q_{lm}$ over every
neighbour inside the cutoff. `dseams.CUTOFF` and `dseams.MUTUAL` select
the neighbour graph. `dseams.set_num_threads` and `dseams.num_threads`
configure OpenMP threads. `dseams.version` returns the engine version.
`dseams.Frame.set_frame` stores one configuration and `dseams.Frame.count`
returns the number of points.

## Ownership

In the `seams` ctypes wrapper, `System`, `Neighbours`, `Rings`,
`CageSet`, `Fingerprint`, and `Brine` are context managers. `close()`
frees the handle (repeated calls are safe), and `__del__` releases
unclosed handles. Neighbour, ring, cage, fingerprint, and brine handles
keep independent copies and remain valid after `System` is closed.
`Neighbours.copy` fills caller-owned NumPy arrays.

## Errors

```python
from seams import SeamsError, System

with System() as sys:
    try:
        sys.chill_plus()
    except SeamsError as err:
        print(err.status, err.message)
        assert err.status == 5
```

An empty frame raises status `5` (`SEAMS_ERR_STATE`). A box that is
neither 3 nor 9 numbers raises `ValueError` before calling the C library.

## Analyses

| Call | Result |
| --- | --- |
| `chill`, `chill_plus` | `ChillResult.labels`, `counts` (length 9) |
| `ql`, `steinhardt` | `OrderResult.ql`, `qlbar` |
| `q4q6`, `q6bar` | `Q4Q6`, `Q6Bar` |
| `neighbours`, `knn` | `Neighbours`, then `NeighbourCSR` |
| `rings` | histogram, labels, vertices, fraction |
| `cages`, `signature` | `CageCounts`, `Signature` |
| `symmetry` | 0 Ih, 1 Ic, 2 stacking, -1 other |
| `rdf`, `running_cn`, `coordination` | `Radial`, `RunningCN`, float |
| `hbonds`, `density` | count, `Density` |
| `fingerprint` | key text and ring census |
| `contact_pairs`, `domains`, `prisms`, `monolayer`, `site_table` | dataclasses or arrays |
| `brine` | fronts, rejection, profile, pockets |

`CUTOFF` is `0`, `MUTUAL` is `1`, and `SEEDED` is `2`. `PUBLISHED` is
`0` and `OVITO` is `1`. `q4q6`, `q6bar`, `ql`, and `steinhardt` keep
scratch buffers per calling thread, so separate frames can run in
separate threads. `set_num_threads` is process-global.

## Adapters and utilities

`dseams.from_ase(atoms)` rotates an ASE cell into the LAMMPS frame and
returns `(positions, box)` where `box` is the nine-number cell vector.
`dseams.from_mdanalysis(universe, select="all")` reads
`universe.dimensions` and atom positions into the same format.

`cages`, `signature_cages`, and `symmetry` use `libseams_lang` (bundled
in the wheel, or pointed to by `SEAMS_LANG_LIB` for custom builds).
`chill_plus`, `q4q6`, `neighbours`, `rings`, and `dseams.Frame` call
`libseams`.

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### `CUTOFF`

- **Purpose**: Graph id of the cutoff neighbour list (`0`).
- **Parameters**: none
- **Returns**: Integer constant.
- **Units**: dimensionless

### `MUTUAL`

- **Purpose**: Graph id of the mutual k-nearest neighbour list (`1`).
- **Parameters**: none
- **Returns**: Integer constant.
- **Units**: dimensionless

### `OVITO`

- **Purpose**: CHILL+ rule (`1`) that averages $q_{lm}$ over every neighbour inside the cutoff.
- **Parameters**: none
- **Returns**: Integer constant.
- **Units**: dimensionless

### `PUBLISHED`

- **Purpose**: CHILL+ rule (`0`) that averages $q_{lm}$ over the four closest neighbours.
- **Parameters**: none
- **Returns**: Integer constant.
- **Units**: dimensionless

### `dseams.Frame`

- **Purpose**: Single-snapshot handle holding coordinates and box dimensions for in-place `dseams` classifiers.
- **Parameters**: none
- **Returns**: A `dseams.Frame` instance.
- **Ownership**: Owns its underlying C `seams_system` and frees it on destruction.
- **Errors**: Methods raise `RuntimeError` with the engine error string when a C call fails.

### `dseams.Frame.set_frame`

- **Purpose**: Store positions and simulation cell on a `dseams.Frame`.
- **Parameters**: `positions` `(n, 3)`, `box` (3 edge lengths or 9 LAMMPS numbers).
- **Returns**: `None`.
- **Ownership**: Copies coordinates into the frame's `seams_system`.
- **Errors**: Raises `ValueError` if `box` is not 3 or 9 numbers; `RuntimeError` if the C call fails.
- **Units**: coordinates in angstrom

### `dseams.Frame.count`

- **Purpose**: Number of points stored on a `dseams.Frame`.
- **Parameters**: `self`
- **Returns**: Integer atom count.
- **Units**: dimensionless count

### `dseams.Frame.chill_plus`

- **Purpose**: CHILL+ labels and nine population counts for the frame.
- **Parameters**: `self`, `cutoff`, `graph`, `rule`
- **Returns**: `(labels, counts)` NumPy arrays owned by the caller.
- **Errors**: `RuntimeError` on engine failure.
- **Units**: `cutoff` in angstrom

### `dseams.Frame.q4q6`

- **Purpose**: Local Steinhardt $q_4$ and $q_6$ on the cutoff shell.
- **Parameters**: `self`, `cutoff`
- **Returns**: `(q4, q6)` NumPy arrays owned by the caller.
- **Errors**: `RuntimeError` on engine failure.
- **Units**: `cutoff` in angstrom; `q4` and `q6` are dimensionless

### `dseams.Frame.neighbours`

- **Purpose**: CSR neighbour list of a `dseams.Frame`, omitting self entries.
- **Parameters**: `self`, `cutoff`, `graph`
- **Returns**: `(offsets, columns)` NumPy arrays owned by the caller.
- **Errors**: `RuntimeError` on engine failure.
- **Units**: `cutoff` in angstrom

### `dseams.Frame.rings`

- **Purpose**: Primitive ring sizes, symmetry labels, size histogram, and six-ring fraction.
- **Parameters**: `self`, `cutoff`, `max_depth`
- **Returns**: `(sizes, symmetry, histogram, fraction)`.
- **Errors**: `RuntimeError` on engine failure.
- **Units**: `cutoff` in angstrom

### `cages`

- **Purpose**: Score hexagonal and double-diamond cage membership.
- **Parameters**: `positions`, `box`, `cutoff=3.5`, `graph=CUTOFF`, `k=4`, `ring_adjacent=False`
- **Returns**: `(per_atom, (n_hex, n_cubic, n_water, n_both))`.
- **Errors**: `ValueError` on invalid `box`; `RuntimeError` on engine failure.
- **Units**: coordinates and `cutoff` in angstrom; counts are dimensionless

### `chill_plus`

- **Purpose**: Label each atom with CHILL+.
- **Parameters**: `positions`, `box`, `rule=PUBLISHED`, `graph=CUTOFF`, `cutoff=3.5`
- **Returns**: `(labels, counts)` NumPy arrays owned by the caller.
- **Errors**: `ValueError` on invalid `box`; `RuntimeError` on engine failure.
- **Units**: coordinates and `cutoff` in angstrom

### `from_ase`

- **Purpose**: Convert an ASE `Atoms` object into positions and a 9-element LAMMPS cell vector.
- **Parameters**: `atoms`
- **Returns**: `(positions, box)` NumPy arrays.
- **Errors**: `ValueError` on degenerate cell vectors.
- **Units**: positions and box lengths in angstrom

### `from_mdanalysis`

- **Purpose**: Extract positions and a 9-element LAMMPS cell vector from an MDAnalysis universe or atom group.
- **Parameters**: `universe`, `select="all"`
- **Returns**: `(positions, box)` NumPy arrays.
- **Errors**: `ValueError` if `dimensions` has fewer than 6 values or zero volume.
- **Units**: positions in angstrom; cell angles in degrees

### `neighbours`

- **Purpose**: Build a cutoff or mutual k-nearest neighbour list as CSR arrays.
- **Parameters**: `positions`, `box`, `cutoff=3.5`, `graph=CUTOFF`
- **Returns**: `(offsets, columns)` NumPy arrays.
- **Errors**: `ValueError` on invalid `box`; `RuntimeError` on engine failure.
- **Units**: coordinates and `cutoff` in angstrom

### `num_threads`

- **Purpose**: Query the OpenMP thread count used by `dseams`.
- **Parameters**: none
- **Returns**: Integer thread count.
- **Units**: dimensionless thread count

### `set_num_threads`

- **Purpose**: Set the process-global OpenMP thread count for `dseams` (`0` restores the library default).
- **Parameters**: `n` (`int`)
- **Returns**: `None`.
- **Units**: dimensionless thread count

### `q4q6`

- **Purpose**: Local Steinhardt $q_4$ and $q_6$ on the cutoff shell.
- **Parameters**: `positions`, `box`, `cutoff=3.5`
- **Returns**: `(q4, q6)` NumPy arrays.
- **Errors**: `ValueError` on invalid `box`; `RuntimeError` on engine failure.
- **Units**: coordinates and `cutoff` in angstrom; `q4` and `q6` are dimensionless

### `rings`

- **Purpose**: Compute primitive-ring sizes, symmetry, histogram, and six-ring fraction.
- **Parameters**: `positions`, `box`, `cutoff=3.5`, `max_depth=6`
- **Returns**: `(sizes, symmetry, histogram, fraction)`.
- **Errors**: `ValueError` on invalid `box`; `RuntimeError` on engine failure.
- **Units**: coordinates and `cutoff` in angstrom; `max_depth` is dimensionless

### `signature_cages`

- **Purpose**: Count closed cages matching a ring-census specification string (such as `"512"`).
- **Parameters**: `positions`, `box`, `spec`, `cutoff=3.5`, `graph=CUTOFF`
- **Returns**: `(n_cages, n_atoms)` tuple of integers.
- **Errors**: `ValueError` on invalid `box`; `RuntimeError` on engine failure.
- **Units**: coordinates and `cutoff` in angstrom

### `symmetry`

- **Purpose**: Shell classification code (`0` Ih, `1` Ic, `2` stacking-disordered, `-1` other).
- **Parameters**: `positions`, `box`, `cutoff=3.5`
- **Returns**: Integer NumPy array of length `n_atoms`.
- **Errors**: `ValueError` on invalid `box`; `RuntimeError` on engine failure.
- **Units**: coordinates and `cutoff` in angstrom

### `version`

- **Purpose**: Return the underlying engine version string.
- **Parameters**: none
- **Returns**: Version string (`str`).

### `userframe.Frame`

- **Purpose**: File-backed frame handle backed by `libseams_lang`.
- **Parameters**: `ptr`, `symbols=None`, `cutoff=3.5`
- **Ownership**: Owns the `seams_frame` pointer; `close()` or deletion frees it.
- **Errors**: Raises `RuntimeError` when `ptr` is null.

### `userframe.read`

- **Purpose**: Read one trajectory frame with `seams_frame_read`.
- **Parameters**: `filename`, `frame=1`, `atom_type=None`, `all_atoms=False`
- **Returns**: A `userframe.Frame`.
- **Errors**: `RuntimeError` on I/O or parser failure.

### `userframe.from_ase`

- **Purpose**: Create a `userframe.Frame` from an ASE `Atoms` selection.
- **Parameters**: `atoms`, `select="O"`, `cutoff=None`, `bonded="auto"`
- **Returns**: A `userframe.Frame`.
- **Errors**: `TypeError` if `atoms` has no positions; `ValueError` if `select` matches nothing.

### `userframe.Frame.n_atoms`

- **Purpose**: Number of atoms stored on the frame.
- **Returns**: Integer atom count.

### `userframe.Frame.positions`

- **Purpose**: Copy of the atomic coordinates stored on the frame.
- **Returns**: `(n, 3)` NumPy float64 array in angstrom.

### `userframe.Frame.box`

- **Purpose**: Simulation cell stored on the frame.
- **Returns**: Three edge lengths, or nine LAMMPS numbers when tilts are non-zero, in angstrom.

### `userframe.Frame.chill_plus`

- **Purpose**: CHILL+ histogram for this frame.
- **Parameters**: `self`, `graph="cutoff"`
- **Returns**: `userframe.IceCounts` mapping label name to population count.
- **Errors**: `RuntimeError` on engine failure.

### `userframe.Frame.cages`

- **Purpose**: Hexagonal and double-diamond cage membership.
- **Parameters**: `self`, `seeded=True`, `k=4`, `ring_adjacent=False`
- **Returns**: `userframe.CageScore`.
- **Errors**: `RuntimeError` on engine failure.

### `userframe.Frame.cages_by_signature`

- **Purpose**: Enumerate closed polyhedra matching a ring-census signature.
- **Parameters**: `self`, `signature`, `max_ring_size=None`, `edge=None`
- **Returns**: List of dicts (`{"signature": ..., "vertices": [...]}`) with `.n_atoms` giving the distinct vertex count.
- **Errors**: `RuntimeError` on engine failure.

### `userframe.Frame.steinhardt`

- **Purpose**: Local $q_l$ and Lechner-Dellago averaged $\bar{q}_l$.
- **Parameters**: `self`, `order_l=6`
- **Returns**: Dictionary with `"ql"` and `"ql_bar"` lists.
- **Errors**: `RuntimeError` on engine failure.

### `userframe.Frame.ring_census`

- **Purpose**: Primitive-ring counts grouped by ring size.
- **Parameters**: `self`, `max_size=8`, `graph=None`
- **Returns**: Dictionary `{size: count}` (zero counts omitted).
- **Errors**: `RuntimeError` on engine failure.

### `userframe.Frame.rdf_2d`

- **Purpose**: In-plane 2D radial distribution function $g(r)$.
- **Parameters**: `self`, `output_dir`, `cutoff=12.0`, `binwidth=0.05`
- **Returns**: `(radii, values)` lists; also writes output under `output_dir`.
- **Errors**: `RuntimeError` on engine failure.

### `userframe.Frame.monolayer_rings`

- **Purpose**: Quasi-2D monolayer ring counts and areal coverage.
- **Parameters**: `self`, `output_dir`, `sheet_area`, `max_depth=4`
- **Returns**: Dictionary `{size: {"count": n, "coverage_xy": area}}`; writes `topoMonolayer` under `output_dir`.
- **Errors**: `RuntimeError` on engine failure.

### `userframe.Frame.find_prisms`

- **Purpose**: Quasi-1D prism census.
- **Parameters**: `self`, `output_dir="output/"`, `max_depth=6`, `shape_matching=False`
- **Returns**: List of `{"size": ..., "perfect": ..., "imperfect": ...}` rows.
- **Errors**: `RuntimeError` on engine failure.

### `userframe.Frame.to_ase`

- **Purpose**: Convert this frame to an ASE `Atoms` object with `ice_type`, `hc`, and `ddc` arrays when computed.
- **Returns**: `ase.Atoms`.
- **Errors**: `ImportError` when ASE is not installed.

### `userframe.Frame.to_solvis` / `userframe.to_solvis`

- **Purpose**: Wrap the frame as a `solvis.System` via `to_ase`.
- **Parameters**: `expand_box=True`
- **Returns**: `solvis.system.System`.
- **Errors**: `ImportError` when `solvis-tools` is not installed.

### `userframe.Frame.close`

- **Purpose**: Free the underlying `seams_frame` handle immediately (safe to call multiple times).
- **Returns**: `None`.

### `userframe.IceCounts`

- **Purpose**: Dictionary subclass holding CHILL+ population counts; missing labels return `0`.

### `userframe.CageScore`

- **Purpose**: Per-atom hexagonal-cage (`hc`) and double-diamond-cage (`ddc`) boolean flags, with dimensionless count properties `n_ih`, `n_ic`, `n_water`, and `n_both`.
<!-- API_SYMBOLS_END -->
