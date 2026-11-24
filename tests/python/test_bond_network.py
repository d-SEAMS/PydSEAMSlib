"""Check graph components and sampled residence counts."""
import numpy as np
import pytest

from pydseams.bond_network import components, residence_counts


def test_components_keep_isolated_sites_and_ignore_edge_direction():
    labels, sizes = components(5, [(2, 1), (0, 1), (1, 2), (4, 4)])
    np.testing.assert_array_equal(labels, [0, 0, 0, 1, 2])
    np.testing.assert_array_equal(sizes, [3, 1, 1])
    assert not components(0, [])[0].size
    with pytest.raises(ValueError, match="site index"):
        components(2, [(0, 2)])


def test_residence_distinguishes_leave_and_return():
    in_region = np.array([[True, False], [False, True], [True, True], [True, True]])
    np.testing.assert_array_equal(
        residence_counts(in_region, [0, 1, 2, 3]),
        [[0, 6, 6, 6], [1, 4, 3, 3], [2, 2, 1, 2], [3, 1, 0, 1]],
    )
    with pytest.raises(ValueError, match="lag"):
        residence_counts(in_region, [4])
