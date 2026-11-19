# dseams tutorial

The `dseams` module exposes a zero-copy NumPy array interface over the
d-SEAMS C API (`seams.h`), alongside `dseams.userframe` for trajectory
file workflows.

## Classify an eight-atom cubic ice cell

`dseams` reads C-contiguous NumPy arrays directly in place. On the
eight-atom diamond cubic unit cell (nearest-neighbour spacing 2.75 Å),
every oxygen is cubic (label `0`) and mean $q_6$ is `0.628539`:

```python
import math
import numpy as np
import dseams

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
pos = np.array([[x * edge, y * edge, z * edge] for x, y, z in sites])
box = (edge, edge, edge)

frame = dseams.Frame()
frame.set_frame(pos, box)
labels, counts = frame.chill_plus(3.5, dseams.CUTOFF, dseams.PUBLISHED)
q4, q6 = frame.q4q6(3.5)
offsets, cols = frame.neighbours(3.5, dseams.CUTOFF)

assert frame.count() == 8
assert int(counts[0]) == 8
assert abs(float(q6.mean()) - 0.628539) < 1e-5
assert len(offsets) == 9 and len(cols) == 32
print("cubic:", int(counts[0]), "mean q6:", round(float(q6.mean()), 6))
```

## Rings, cages, and shell symmetry on a supercell

Topological cage and primitive-ring analyses require periodic boxes
wider than twice the neighbour cutoff (`> 7.0` Å). On a 3x3x3 cubic
diamond supercell (216 atoms), every atom belongs to a double-diamond
cage (`n_cubic == 216`) and every primitive ring is a six-ring:

```python
import math
import numpy as np
import dseams

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
pos = np.array(
    [
        [(i + x) * edge, (j + y) * edge, (k + z) * edge]
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for x, y, z in sites
    ]
)
box = (3 * edge, 3 * edge, 3 * edge)

sizes, symmetry, hist, fraction = dseams.rings(pos, box)
per_atom, (n_hex, n_cubic, n_water, n_both) = dseams.cages(pos, box)
n_cages, n_atoms = dseams.signature_cages(pos, box, spec="512")
shell = dseams.symmetry(pos, box)

assert int(hist[6]) == 432 and fraction == 1.0
assert (n_hex, n_cubic, n_water, n_both) == (0, 216, 0, 0)
assert (n_cages, n_atoms) == (0, 0)
assert np.all(shell == 1)
print("six-rings:", int(hist[6]), "DDC atoms:", n_cubic)
```

## Convert ASE or MDAnalysis frames

`dseams.from_ase(atoms)` rotates any nonsingular periodic ASE cell into
LAMMPS restricted-triclinic coordinates, and `dseams.from_mdanalysis`
extracts positions and cell parameters from an MDAnalysis universe or
atom group:

```python
from ase import Atoms
import dseams

atoms = Atoms("O" * len(pos), positions=pos, cell=list(box), pbc=True)
apos, abox = dseams.from_ase(atoms)
labels, counts = dseams.chill_plus(apos, abox, rule=dseams.OVITO)
assert len(labels) == len(atoms)
assert int(counts[0]) == 216
```

You can also run the built-in cubic-cell self-check from the command
line:

```bash
python -m dseams
```
