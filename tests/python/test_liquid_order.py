import numpy as np
import pytest

from pydseams import Frame


def test_tetrahedral_q_in_skewed_periodic_cell():
    lattice = np.array([[12.0, 0.0, 0.0], [2.0, 11.0, 0.0], [1.0, 3.0, 10.0]])
    center = np.array([0.4, 0.5, 0.6])
    directions = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(
        3
    )
    vertices = (
        center
        + 2.75 * directions
        + np.array([lattice[1], -lattice[2], lattice[0], -lattice[1]])
    )
    frame = Frame.from_arrays(
        np.vstack((center, vertices)), [15.0, 14.0, 10.0, 2.0, 1.0, 3.0]
    )
    q = frame.tetrahedral_q()
    assert q[0] == pytest.approx(1.0)
    assert np.isnan(q[1:]).all()
    assert np.isfinite(frame.tetrahedral_q(cutoff=None)).all()


def test_lsi_includes_first_outside_neighbor():
    frame = Frame.from_arrays(
        [[x, 1.0, 1.0] for x in (0.0, 1.0, 2.0, 3.0, 4.4)], [20.0] * 3
    )
    assert frame.local_structure_index()[0] == pytest.approx(np.var([1.0, 1.0, 1.4]))
    assert np.isnan(frame.local_structure_index(cutoff=0.9)).all()
    assert np.isnan(frame.local_structure_index(cutoff=6.0)).all()
