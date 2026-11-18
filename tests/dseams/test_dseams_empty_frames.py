# Copyright (c) d-SEAMS developers
"""Empty frames on dseams, userframe, and the ctypes seams package.

Zero atoms are rejected, or they return an empty array. A dump whose
oxygen guess keeps nothing is a frame with no molecules. Three atoms
farther apart than the cutoff have an empty neighbour list.
"""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms

import dseams
from dseams import frames, userframe
from seams import CUTOFF, ERR_STATE, SeamsError, System

BOX = (30.0, 30.0, 30.0)
ISOLATED = np.array(
    [[1.0, 1.0, 1.0], [12.0, 1.0, 1.0], [1.0, 12.0, 1.0]], dtype=np.float64
)
EMPTY = np.zeros((0, 3), dtype=np.float64)
EMPTY_DUMP = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
0
ITEM: BOX BOUNDS pp pp pp
0 30
0 30
0 30
ITEM: ATOMS id mol type x y z
"""
NO_OXYGEN = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0 30
0 30
0 30
ITEM: ATOMS id mol type x y z
1 0 3 1.0 1.0 1.0
2 0 3 12.0 1.0 1.0
3 0 3 1.0 12.0 1.0
"""


def _i32(array, shape):
    assert array.shape == shape
    assert array.dtype == np.int32


def _f64(array, shape):
    assert array.shape == shape
    assert array.dtype == np.float64


def test_nanobind_rejects_zero_atoms_and_an_unset_frame():
    frame = dseams.Frame()
    box = np.asarray(BOX, dtype=np.float64)
    with pytest.raises(RuntimeError, match="system has no frame"):
        frame.chill_plus()
    with pytest.raises(RuntimeError, match="system has no frame"):
        frame.q4q6()
    with pytest.raises(RuntimeError, match="system has no frame"):
        frame.neighbours()
    with pytest.raises(RuntimeError, match="system has no frame"):
        frame.rings()
    assert frame.count() == 0
    with pytest.raises(RuntimeError, match="frame has no atoms"):
        frame.set_frame(EMPTY, box)
    with pytest.raises(RuntimeError, match="frame has no atoms"):
        dseams.chill_plus(EMPTY, BOX)
    with pytest.raises(RuntimeError, match="frame has no atoms"):
        dseams.q4q6(EMPTY, BOX)
    with pytest.raises(RuntimeError, match="frame has no atoms"):
        dseams.neighbours(EMPTY, BOX)
    with pytest.raises(RuntimeError, match="frame has no atoms"):
        dseams.rings(EMPTY, BOX)
    with pytest.raises(RuntimeError, match="positive atom count"):
        frames.rings(EMPTY, BOX)
    with pytest.raises(RuntimeError, match="positive atom count"):
        frames.cages(EMPTY, BOX)
    with pytest.raises(RuntimeError, match="positive atom count"):
        frames.signature_cages(EMPTY, BOX)
    with pytest.raises(RuntimeError, match="positive atom count"):
        frames.symmetry(EMPTY, BOX)


def test_nanobind_isolated_atoms_are_water_with_no_neighbours():
    frame = dseams.Frame()
    frame.set_frame(ISOLATED, np.asarray(BOX, dtype=np.float64))
    assert frame.count() == 3
    labels, counts = frame.chill_plus()
    _i32(labels, (3,))
    _i32(counts, (9,))
    assert np.all(labels == 2)
    assert int(counts[2]) == 3 and int(counts.sum()) == 3
    q4, q6 = frame.q4q6()
    _f64(q4, (3,))
    _f64(q6, (3,))
    assert np.all(q4 == 0.0) and np.all(q6 == 0.0)
    offsets, columns = frame.neighbours()
    _i32(offsets, (4,))
    _i32(columns, (0,))
    assert np.all(offsets == 0)
    sizes, sym, hist, fraction = frame.rings()
    _i32(sizes, (0,))
    _i32(sym, (0,))
    _i32(hist, (13,))
    assert np.all(hist == 0) and fraction == 0.0
    per, populations = frames.cages(ISOLATED, BOX, graph="cutoff")
    _i32(per, (3,))
    assert np.all(per == 0)
    assert populations == (0, 0, 3, 0)
    assert frames.signature_cages(ISOLATED, BOX) == (0, 0)
    symmetry = frames.symmetry(ISOLATED, BOX)
    _i32(symmetry, (3,))
    assert np.all(symmetry == -1)
    free_labels, free_counts = dseams.chill_plus(ISOLATED, BOX)
    assert np.array_equal(free_labels, labels)
    assert np.array_equal(free_counts, counts)


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
            # The engine still reports pockets when every ice flag is clear.
            assert pockets.n_pockets == 2
            assert pockets.reservoir_waters == 1


def _empty_userframe(path, text, **kwargs):
    path.write_text(text, encoding="utf-8")
    frame = userframe.read(path, **kwargs)
    assert frame.n_atoms == 0
    positions = frame.positions
    _f64(positions, (0, 3))
    assert frame.chill_plus().water == 0
    assert frame.ring_census() == {}
    order = frame.steinhardt()
    assert order["ql"] == [] and order["ql_bar"] == []
    cages = frame.cages(seeded=False)
    assert cages.n_ih == 0 and cages.n_ic == 0 and cages.n_water == 0
    assert frame.cages_by_signature("512") == []
    centres, values = frame.rdf_2d(path.parent / "rdf2d", cutoff=6.0, binwidth=0.5)
    assert centres == [0.25 + 0.5 * i for i in range(12)]
    assert values == [0.0] * 12
    rows = frame.find_prisms(path.parent / "prisms")
    assert rows == [
        {"size": size, "perfect": 0, "imperfect": 0} for size in range(3, 7)
    ]
    monolayer = frame.monolayer_rings(path.parent / "mono", 100.0)
    assert monolayer[3]["count"] == 0 and monolayer[4]["count"] == 0
    assert monolayer[3]["coverage_xy"] == 0.0
    atoms = frame.to_ase()
    assert len(atoms) == 0
    assert atoms.positions.shape == (0, 3)


def test_userframe_zero_atoms_and_zero_oxygen(tmp_path):
    _empty_userframe(tmp_path / "empty.lammpstrj", EMPTY_DUMP)
    _empty_userframe(tmp_path / "no-oxygen.lammpstrj", NO_OXYGEN)
    atoms = Atoms(cell=BOX, pbc=True)
    with pytest.raises(ValueError, match="no atoms matched"):
        userframe.from_ase(atoms, select=None)


def test_userframe_isolated_atoms(tmp_path):
    text = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0 30
0 30
0 30
ITEM: ATOMS id mol type x y z
1 1 1 1.0 1.0 1.0
2 1 1 12.0 1.0 1.0
3 1 1 1.0 12.0 1.0
"""
    path = tmp_path / "far.lammpstrj"
    path.write_text(text, encoding="utf-8")
    frame = userframe.read(path, atom_type=1)
    assert frame.n_atoms == 3
    _f64(frame.positions, (3, 3))
    assert np.all(frame.positions == ISOLATED)
    assert frame.chill_plus().water == 3
    assert frame.ring_census() == {}
    order = frame.steinhardt()
    assert order["ql"] == [0.0, 0.0, 0.0]
    assert order["ql_bar"] == [0.0, 0.0, 0.0]
    assert frame.cages(seeded=False).n_water == 3
    hydrogens = Atoms("H3", positions=ISOLATED, cell=BOX, pbc=True)
    adopted = userframe.from_ase(hydrogens, select="H", cutoff=3.5)
    assert adopted.n_atoms == 3
    assert adopted.chill_plus().water == 3
    assert adopted.ring_census() == {}
