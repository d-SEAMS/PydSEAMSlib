"""d-SEAMS ice states as an OVITO Python modifier.

Register :func:`ice_states` with ``pipeline.modifiers.append(ice_states)``
(OVITO 3.x ``PythonScriptModifier`` function form). It adds the particle
property ``Ice state`` (0 water, 1 cubic, 2 hexagonal, 3 mixed) computed
by the seeded cage assignment on the oxygens of the selected type, so the
Color coding modifier paints ice growth frame by frame. OVITO itself is
never imported here: the function only reads the ``DataCollection`` it is
handed.
"""

from __future__ import annotations

import numpy as np

from pydseams.features import IceFeaturizer
from pydseams.frame import Frame

PROPERTY = "Ice state"


def _cell_box(cell):
    m = np.asarray(cell, dtype=float)
    h = m[:, :3]
    origin = np.array(
        [float(x) for x in m[:, 3]] if m.shape[1] > 3 else [0.0, 0.0, 0.0],
        dtype=float,
    )
    # OVITO SimulationCell columns are a, b, c. Encode the same LAMMPS
    # dump bound spans as aseio._dump_box so dumpBoundsToH recovers H.
    lx = float(h[0, 0])
    xy = float(h[0, 1])
    xz = float(h[0, 2])
    ly = float(h[1, 1])
    yz = float(h[1, 2])
    lz = float(h[2, 2])
    if max(abs(xy), abs(xz), abs(yz)) < 1e-9:
        return [lx, ly, lz], origin.tolist()
    xmin = min(0.0, xy, xz, xy + xz)
    xmax = max(0.0, xy, xz, xy + xz)
    ymin = min(0.0, yz)
    ymax = max(0.0, yz)
    box = [lx + xmax - xmin, ly + ymax - ymin, lz, xy, xz, yz]
    box_low = (origin + np.array([xmin, ymin, 0.0])).tolist()
    return box, box_low


def ice_states(frame, data, oxygen_type=None, cutoff=3.5, k=4, ring_adjacent=True):
    """Modifier function: adds the ``Ice state`` particle property.

    ``oxygen_type`` selects the particles by their ``Particle Type`` id;
    ``None`` takes every particle. Particles outside the selection get
    ``-1``.
    """
    particles = data.particles
    positions = np.asarray(particles.positions, dtype=float)
    if oxygen_type is None:
        mask = np.ones(len(positions), dtype=bool)
    else:
        mask = np.asarray(particles.particle_types) == int(oxygen_type)
    box, box_low = _cell_box(data.cell)
    fr = Frame.from_arrays(positions[mask], box, cutoff=cutoff, box_low=box_low)
    _, states = IceFeaturizer(
        fr, k=k, ring_adjacent=ring_adjacent, chill=False
    ).frame_features()
    out = np.full(len(positions), -1, dtype=np.int32)
    out[mask] = np.asarray(states, dtype=np.int32)
    data.particles_.create_property(PROPERTY, data=out)
