# Copyright (c) d-SEAMS developers
"""Empty-frame checks for the source-only ctypes seams package."""

from __future__ import annotations

import numpy as np
import pytest

from seams import CUTOFF, ERR_STATE, SeamsError, System

BOX = (30.0, 30.0, 30.0)
ISOLATED = np.array(
    [[1.0, 1.0, 1.0], [12.0, 1.0, 1.0], [1.0, 12.0, 1.0]], dtype=np.float64
)
EMPTY = np.zeros((0, 3), dtype=np.float64)


def _i32(array, shape):
    assert array.shape == shape
    assert array.dtype == np.int32


def _f64(array, shape):
    assert array.shape == shape
    assert array.dtype == np.float64


def test_seams_rejects_empty_positions_and_an_unset_system():
    with System() as system:
        with pytest.raises(ValueError, match="positions are empty"):
            system.set_frame(EMPTY, BOX)
        with pytest.raises(ValueError, match="positions are empty"):
            system.set_ortho(EMPTY, BOX)
        assert system.count() == 0
        calls = [
            lambda: system.neighbours(),
            lambda: system.knn(),
            lambda: system.chill(),
            lambda: system.chill_plus(),
            lambda: system.ql(6),
            lambda: system.steinhardt(6),
            lambda: system.q4q6(),
            lambda: system.q6bar(),
            lambda: system.rings(),
            lambda: system.cages(),
            lambda: system.signature(),
            lambda: system.symmetry(),
            lambda: system.rdf(),
            lambda: system.running_cn(),
            lambda: system.coordination(),
            lambda: system.hbonds(),
            lambda: system.density(),
            lambda: system.fingerprint(),
            lambda: system.prisms(),
            lambda: system.monolayer(),
            lambda: system.contact_pairs("8=cationHead,9=anion"),
            lambda: system.domains("8=cationHead,9=anion"),
            lambda: system.site_table("1=waterO"),
            lambda: system.brine(np.zeros(0, dtype=np.int32)),
        ]
        for call in calls:
            with pytest.raises(SeamsError, match="no frame") as caught:
                call()
            assert caught.value.status == ERR_STATE


def test_seams_isolated_atoms_have_an_empty_cutoff_graph():
    with System() as system:
        system.set_ortho(ISOLATED, BOX)
        assert system.count() == 3
        with system.neighbours(cutoff=3.5, graph=CUTOFF) as neighbours:
            csr = neighbours.copy()
        _i32(csr.offsets, (4,))
        _i32(csr.columns, (0,))
        assert np.all(csr.offsets == 0)
        with system.knn(k=4, cutoff=3.5) as neighbours:
            knn = neighbours.copy()
        _i32(knn.offsets, (4,))
        assert knn.columns.dtype == np.int32
        rows = [
            set(knn.columns[int(knn.offsets[i]) : int(knn.offsets[i + 1])])
            for i in range(3)
        ]
        assert rows == [{1, 2}, {0, 2}, {0, 1}]

        published = system.chill_plus(cutoff=3.5)
        _i32(published.labels, (3,))
        _i32(published.counts, (9,))
        assert np.all(published.labels == 2)
        assert int(published.counts[2]) == 3 and int(published.counts.sum()) == 3
        original = system.chill(cutoff=3.5)
        _i32(original.labels, (3,))
        assert int(original.counts.sum()) == 3

        with system.rings(cutoff=3.5) as rings:
            assert rings.count() == 0
            hist = rings.histogram()
            _i32(hist, (13,))
            assert np.all(hist == 0)
            labels = rings.labels()
            _i32(labels.sizes, (0,))
            _i32(labels.symmetry, (0,))
            assert rings.fraction(6) == 0.0

        order = system.q4q6(cutoff=3.5)
        _f64(order.q4, (3,))
        _f64(order.q6, (3,))
        assert np.all(order.q4 == 0.0) and np.all(order.q6 == 0.0)
        local = system.ql(6, cutoff=3.5)
        _f64(local.ql, (3,))
        assert local.qlbar is None and np.all(local.ql == 0.0)
        averaged = system.steinhardt(6, cutoff=3.5, average=True)
        _f64(averaged.ql, (3,))
        _f64(averaged.qlbar, (3,))
        assert np.all(averaged.ql == 0.0) and np.all(averaged.qlbar == 0.0)
        bar = system.q6bar(cutoff=3.5)
        assert bar.n_ice == 0 and bar.mean == 0.0 and bar.median == 0.0
        _f64(bar.values, (3,))
        assert np.all(bar.values == 0.0)

        with system.cages(cutoff=3.5, graph=CUTOFF) as cages:
            counts = cages.counts()
            assert (counts.n_hex, counts.n_cubic, counts.n_water, counts.n_both) == (
                0,
                0,
                3,
                0,
            )
            cage_labels = cages.labels(3)
            _i32(cage_labels, (3,))
            assert np.all(cage_labels == 0)
        signature = system.signature("512", cutoff=3.5)
        assert signature.n_cages == 0 and signature.n_atoms == 0
        symmetry = system.symmetry()
        _i32(symmetry, (3,))
        assert np.all(symmetry == -1)

        radial = system.rdf(rmax=6.0, bins=4)
        _f64(radial.r, (4,))
        _f64(radial.g, (4,))
        _i32(radial.count, (4,))
        assert np.all(radial.g == 0.0) and np.all(radial.count == 0)
        assert radial.volume == 27000.0
        running = system.running_cn(rmax=6.0, bins=4)
        _f64(running.cn, (4,))
        assert np.all(running.cn == 0.0)
        assert system.coordination(rmax=3.5) == 0.0
        assert system.hbonds() == 0
        density = system.density(bins=4)
        _f64(density.z, (4,))
        _f64(density.rho, (4,))
        assert density.rho[0] == 3.0 / (30.0 * 30.0 * (30.0 / 4.0))
        assert np.all(density.rho[1:] == 0.0)
        pairs = system.contact_pairs("8=cationHead,9=anion")
        assert pairs.n_pairs == 0 and pairs.n_cation == 0 and pairs.n_anion == 0
        domains = system.domains("8=cationHead,9=anion")
        assert domains.n_subset == 0 and domains.largest == 0
        assert domains.percolation == 0.0
        sites = system.site_table("1=waterO")
        _i32(sites, (3,))
        monolayer = system.monolayer(cutoff=3.5, max_depth=6, sheet_area=100.0)
        _i32(monolayer, (7,))
        assert np.all(monolayer == 0)
        prisms = system.prisms(cutoff=3.5)
        _i32(prisms.perfect, (7,))
        _i32(prisms.imperfect, (7,))
        assert np.all(prisms.perfect == 0) and np.all(prisms.imperfect == 0)
        with system.fingerprint(cutoff=3.5) as fingerprint:
            assert fingerprint.key()
        with system.brine(np.zeros(3, dtype=np.int32), n_bins=4) as brine:
            fronts = brine.fronts()
            assert fronts.found is False
            profile = brine.profile()
            _f64(profile.z_mid, (4,))
            _f64(profile.molality, (4,))
            _f64(profile.ice_fraction, (4,))
            pockets = brine.pockets()
            # seams-core 2.14.2 counts every water of an all-liquid frame as reservoir.
            assert pockets.n_pockets == 0
            assert pockets.reservoir_waters == 3
