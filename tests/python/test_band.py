"""The pair that leaves a cutoff graph is named on the next image."""

import pytest

from pydseams import Frame
from pydseams.band import band_rows

CELL = [20.0, 20.0, 20.0]


def test_band_rows_names_the_pair_that_leaves():
    side = 1.5
    first = Frame.from_arrays(
        [(0.0, 0.0, 0.0), (side, 0.0, 0.0), (side, side, 0.0), (0.0, side, 0.0)],
        CELL,
        numbers=[1, 1, 1, 1],
        cutoff=2.0,
    )
    second = Frame.from_arrays(
        [(0.0, 0.0, 0.0), (2.5, 0.0, 0.0), (side, side, 0.0), (0.0, side, 0.0)],
        CELL,
        numbers=[1, 1, 1, 1],
        cutoff=2.0,
    )
    rows = band_rows([first, second], max_ring_size=8)
    assert rows[0]["entered"] == []
    assert rows[0]["left"] == []
    assert rows[0]["rings"] == {4: 1}
    assert rows[1]["left"] == [(0, 1)]
    assert (0, 1) not in rows[1]["entered"]
    assert rows[1]["rings"] == {}
    assert rows[1]["n_bonds"] == rows[0]["n_bonds"] - 1


def test_band_rows_records_a_new_same_type_contact():
    early = Frame.from_arrays(
        [(0.0, 0.0, 0.0), (1.6, 0.0, 0.0), (1.6, 3.0, 0.0)],
        CELL,
        numbers=[1, 2, 2],
        cutoff=2.0,
    )
    late = Frame.from_arrays(
        [(0.0, 0.0, 0.0), (2.4, 0.0, 0.0), (2.4, 1.1, 0.0)],
        CELL,
        numbers=[1, 2, 2],
        cutoff=2.0,
    )
    rows = band_rows([early, late])
    assert (0, 1) in rows[1]["left"]
    assert (1, 2) in rows[1]["entered"]
    assert rows[1]["shortest_within_type"][2] == pytest.approx(1.1)


@pytest.mark.parametrize("translation", [0.0, 60.0, -40.0])
def test_band_distance_uses_periodic_images(translation):
    frame = Frame.from_arrays(
        [(0.2, 0.0, 0.0), (19.4 + translation, 0.0, 0.0)],
        CELL, numbers=[1, 1], cutoff=2.0,
    )
    assert band_rows([frame])[0]["shortest_within_type"][1] == pytest.approx(0.8)
    open_frame = Frame(cloud=frame.cloud, periodic=False, bonded="cutoff", cutoff=2.0)
    assert band_rows([open_frame])[0]["shortest_within_type"][1] == pytest.approx(
        abs(19.2 + translation)
    )


def test_band_distance_in_skew_cell_matches_lattice_enumeration():
    from itertools import product
    from math import sqrt

    delta = (1.6, 0.9, 0.0)
    expected = min(
        sqrt((delta[0] - 4 * i - 3 * j) ** 2
             + (delta[1] - 2 * j) ** 2 + (delta[2] - 5 * k) ** 2)
        for i, j, k in product(range(-3, 4), repeat=3)
    )
    frame = Frame.from_arrays(
        [(0.0, 0.0, 0.0), delta], [7.0, 2.0, 5.0, 3.0, 0.0, 0.0],
        numbers=[1, 1], cutoff=2.0,
    )
    assert band_rows([frame])[0]["shortest_within_type"][1] == pytest.approx(expected)


def test_band_rejects_changed_atom_identity():
    first = Frame.from_arrays([(0, 0, 0), (1, 0, 0)], CELL, numbers=[1, 2])
    changed = Frame.from_arrays([(0, 0, 0), (1, 0, 0)], CELL, numbers=[2, 1])
    with pytest.raises(ValueError, match="IDs, types, and order"):
        band_rows([first, changed])


def test_periodic_distance_checks_atom_indices():
    from pydseams import yoda

    frame = Frame.from_arrays([(0, 0, 0), (1, 0, 0)], CELL)
    with pytest.raises(IndexError):
        yoda.periodicDistance(frame.cloud, -1, 0)
    with pytest.raises(IndexError):
        yoda.periodicDistance(frame.cloud, 0, 2)
