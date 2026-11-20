# Copyright (c) d-SEAMS developers
"""Empty frames on the installed dseams and userframe surfaces.

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
