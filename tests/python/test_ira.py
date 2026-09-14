"""IRA/SOFI bindings. Skip the overlay when this build did not link libira."""

import pydseams as ds
import pytest


def test_ira_available_is_bool():
    assert isinstance(ds.yoda.ira_available(), bool)


@pytest.mark.skipif(not ds.yoda.ira_available(), reason="libira not linked")
def test_ira_match_rotated_square():
    ref = [
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
    ]
    tgt = [
        [-1.0, 1.0, 0.0],
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [1.0, 1.0, 0.0],
    ]
    err, rmsd, hausdorff, quat, assignment = ds.yoda.ira_match(ref, tgt)
    assert err == 0
    assert rmsd == pytest.approx(0.0, abs=1e-6)
    assert hausdorff == pytest.approx(0.0, abs=1e-6)
    assert len(quat) == 4
    assert len(assignment) == 4


@pytest.mark.skipif(not ds.yoda.ira_available(), reason="libira not linked")
def test_sofi_square_not_c1():
    pts = [
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
    ]
    err, symbol, n_ops = ds.yoda.sofi_point_group(pts)
    assert err == 0
    assert n_ops >= 2
    assert symbol.strip() != "C1"
