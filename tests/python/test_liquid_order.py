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


def test_c_api_frame_reuses_current_coordinates():
    import dseams

    frame = dseams.Frame()
    xyz = np.array(
        [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4.4, 0, 0]], dtype=float
    )
    box = [1.0, -2.0, 3.0, 20.0, 20.0, 20.0, 5.0, 2.0, 4.0]
    frame.set_frame(xyz, box)
    assert frame.local_structure_index()[0] == pytest.approx(np.var([1.0, 1.0, 1.4]))
    xyz[4, 0] = 4.1
    frame.set_frame(xyz, box)
    assert frame.local_structure_index()[0] == pytest.approx(np.var([1.0, 1.0, 1.1]))
    assert np.isnan(frame.tetrahedral_q()[0])
    assert frame.tetrahedral_q(cutoff=None)[0] == pytest.approx(-3.0)


def test_angular_order_uses_euclidean_image_in_oblique_cell():
    from itertools import product

    cell = np.array([[12.0, 0.0, 0.0], [5.9, 10.0, 0.0], [0.0, 0.0, 12.0]])
    far = np.array([0.49, 0.49, 0.2]) @ cell
    points = np.vstack(([0.0, 0.0, 0.0], far, np.eye(3)))
    candidates = far + np.array(list(product(range(-1, 2), repeat=3))) @ cell
    nearest = candidates[np.argmin(np.linalg.norm(candidates, axis=1))]
    assert np.linalg.norm(nearest) < np.linalg.norm(far)
    unit = np.vstack((nearest / np.linalg.norm(nearest), np.eye(3)))
    cosines = (unit @ unit.T)[np.triu_indices(4, 1)]
    expected = 1 - (3 / 8) * np.sum((cosines + 1 / 3) ** 2)
    frame = Frame.from_arrays(points, [17.9, 10.0, 12.0, 5.9, 0.0, 0.0])
    assert frame.tetrahedral_q(cutoff=None)[0] == pytest.approx(expected)
