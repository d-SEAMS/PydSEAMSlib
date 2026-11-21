import numpy as np
import pytest

from pydseams import Frame


def test_tetrahedral_q_in_skewed_periodic_cell():
    lattice = np.array([[12., 0., 0.], [2., 11., 0.], [1., 3., 10.]])
    center = np.array([.4, .5, .6])
    directions = np.array([[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1]]) / np.sqrt(3)
    vertices = center + 2.75 * directions + np.array([lattice[1], -lattice[2], lattice[0], -lattice[1]])
    frame = Frame.from_arrays(np.vstack((center, vertices)), [15., 14., 10., 2., 1., 3.])
    q = frame.tetrahedral_q()
    assert q[0] == pytest.approx(1.)
    assert np.isnan(q[1:]).all()
    assert np.isfinite(frame.tetrahedral_q(cutoff=None)).all()


def test_lsi_includes_first_outside_neighbor():
    frame = Frame.from_arrays([[x,1.,1.] for x in (0., 1., 2., 3., 4.4)], [20.]*3)
    assert frame.local_structure_index()[0] == pytest.approx(np.var([1., 1., 1.4]))
    assert np.isnan(frame.local_structure_index(cutoff=.9)).all()
    assert np.isnan(frame.local_structure_index(cutoff=6.)).all()
