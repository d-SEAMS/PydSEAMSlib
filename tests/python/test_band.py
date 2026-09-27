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
