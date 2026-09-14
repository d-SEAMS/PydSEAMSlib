import numpy as np
import pytest

from pydseams.adapters import IceStates, ice_states
from pydseams.features import STATE_IC, STATE_WATER

from test_features import _cubic_diamond


class _Particles:
    def __init__(self, pos, types):
        self.positions = pos
        self.particle_types = types
        self.props = {}

    def create_property(self, name, data):
        self.props[name] = np.asarray(data)


class _Data:
    def __init__(self, pos, types, cell):
        self.particles = _Particles(pos, types)
        self.particles_ = self.particles
        self.cell = cell


def test_ovito_modifier_paints_a_cubic_lattice_and_skips_other_types():
    pos, cell = _cubic_diamond(4)
    guest = np.array([[0.7, 0.7, 0.7]])
    allpos = np.vstack([pos, guest])
    types = np.array([1] * len(pos) + [2])
    ovcell = np.column_stack([np.diag(cell), np.zeros(3)])
    data = _Data(allpos, types, ovcell)
    ice_states(0, data, oxygen_type=1, cutoff=3.5)
    out = data.particles.props["Ice state"]
    assert out.shape == (len(pos) + 1,)
    assert out[-1] == -1
    assert (out[:-1] == STATE_IC).mean() > 0.9
    assert (out[:-1] != STATE_WATER).all() or (out[:-1] == STATE_WATER).sum() < len(
        pos
    ) // 10


def test_ovito_cell_box_uses_dump_bound_spans_for_tilt():
    from pydseams.adapters._ovito import _cell_box

    # OVITO columns: a=(10,0,0), b=(2,8,0), c=(1,1,6), origin=0
    h = np.array(
        [
            [10.0, 2.0, 1.0, 0.0],
            [0.0, 8.0, 1.0, 0.0],
            [0.0, 0.0, 6.0, 0.0],
        ]
    )
    box, box_low = _cell_box(h)
    xy, xz, yz = 2.0, 1.0, 1.0
    xmin = min(0.0, xy, xz, xy + xz)
    xmax = max(0.0, xy, xz, xy + xz)
    ymin = min(0.0, yz)
    ymax = max(0.0, yz)
    assert box == [10.0 + xmax - xmin, 8.0 + ymax - ymin, 6.0, xy, xz, yz]
    assert box_low == [xmin, ymin, 0.0]


def test_mdanalysis_ice_states_over_a_two_frame_lattice():
    mda = pytest.importorskip("MDAnalysis")
    pos, cell = _cubic_diamond(4)
    n = len(pos)
    u = mda.Universe.empty(n, trajectory=True)
    u.add_TopologyAttr("names", ["O"] * n)
    u.add_TopologyAttr("types", ["O"] * n)
    from MDAnalysis.coordinates.memory import MemoryReader

    coords = np.array([pos, pos], dtype=np.float32)
    dims = np.array(
        [[cell[0], cell[1], cell[2], 90.0, 90.0, 90.0]] * 2, dtype=np.float32
    )
    u.trajectory = MemoryReader(coords, dimensions=dims)
    an = IceStates(u.atoms, cutoff=3.5).run()
    assert an.results.states.shape == (2, n)
    assert (an.results.states[0] == STATE_IC).mean() > 0.9
    assert an.results.features.shape[0] == 2
    assert "n_ic" in an.results.names or any("ic" in s for s in an.results.names)
