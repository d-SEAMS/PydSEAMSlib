"""Radius and ray-parity occupancy of guests in the two 5^12 cages of sI."""

from pathlib import Path

import numpy as np

from pydseams import yoda
from pydseams.frame import DualOccupancy, Frame


def _si_waters():
    """Oxygen positions and edge of the GenIce sI cell shipped with the tests."""
    path = Path(__file__).resolve().parents[1] / "data" / "genice_sI.lammpstrj"
    lines = path.read_text(encoding="utf-8").splitlines()
    edge = float(lines[5].split()[1]) - float(lines[5].split()[0])
    pos = np.array([[float(v) for v in row.split()[2:5]] for row in lines[9:]])
    return pos, edge


def _periodic_distance(a, b, edge):
    d = np.asarray(a) - np.asarray(b)
    d -= edge * np.round(d / edge)
    return float(np.linalg.norm(d))


def test_guest_occupancy_both_agrees_on_si():
    pos, edge = _si_waters()
    waters = Frame.from_arrays(pos, [edge] * 3, cutoff=3.5)
    found = waters.cages_by_signature("512")
    assert len(found) == 2
    centres = [
        yoda.periodicCentroid(waters.cloud, [int(v) for v in c["vertices"]])
        for c in found
    ]
    # A third guest on the grid point farthest from both 5^12 centres sits in
    # a 5^12 6^2 cage, so neither assignment places it.
    grid = [
        np.array([i, j, k]) * edge / 4.0
        for i in range(4)
        for j in range(4)
        for k in range(4)
    ]
    far = max(grid, key=lambda p: min(_periodic_distance(p, c, edge) for c in centres))
    assert min(_periodic_distance(far, c, edge) for c in centres) > 5.0
    guests = np.vstack([np.asarray(centres), far])
    frame = Frame.from_arrays(
        np.vstack([pos, guests]),
        [edge] * 3,
        numbers=[1] * len(pos) + [2] * 3,
        cutoff=3.5,
    )
    both = frame.guest_occupancy_both("512", guest_types=[2], radius=4.0)
    assert isinstance(both, DualOccupancy)
    assert list(both.radius.cageOfGuest) == [0, 1, -1]
    assert list(both.inside.cageOfGuest) == [0, 1, -1]
    for occ in (both.radius, both.inside):
        assert list(occ.guestsPerCage) == [1, 1]
        assert (occ.occupied, occ.multiply, occ.free) == (2, 0, 1)
