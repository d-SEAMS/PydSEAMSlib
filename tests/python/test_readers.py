"""Real reads of eOn CON, PDB, and GRO.

Each test is skipped only when the loaded engine reports that the reader
was not linked. A build that has the reader must run the assertion.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pydseams import available_readers, has_chemfiles, has_readcon, read

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "readers"
CON = FIXTURES / "tiny_multi_cuh2.con"
PDB = FIXTURES / "water.pdb"
GRO = FIXTURES / "water.gro"


def test_feature_query_matches_reader_flags():
    readers = available_readers()
    assert has_chemfiles() is readers["chemfiles"]
    assert has_readcon() is readers["readcon"]
    assert readers["lammps"] is True


def test_read_con_file():
    if not has_readcon():
        pytest.skip("readcon-core is not linked")
    frame = read(CON)
    assert frame.n_atoms == 4
    assert frame.box[0] == pytest.approx(15.3456)
    assert frame.box[1] == pytest.approx(21.702)
    assert frame.box[2] == pytest.approx(100.0)
    assert len(frame.positions) == 4


def test_read_pdb_file():
    if not has_chemfiles():
        pytest.skip("chemfiles is not linked")
    frame = read(PDB)
    assert frame.n_atoms == 3
    oxygen = frame.positions[0]
    assert oxygen[0] == pytest.approx(5.0)
    assert oxygen[1] == pytest.approx(5.0)
    assert oxygen[2] == pytest.approx(5.0)


def test_read_gro_file():
    if not has_chemfiles():
        pytest.skip("chemfiles is not linked")
    frame = read(GRO)
    assert frame.n_atoms == 3
    oxygen = frame.positions[0]
    assert oxygen[0] == pytest.approx(5.0)
    assert oxygen[1] == pytest.approx(5.0)
    assert oxygen[2] == pytest.approx(5.0)
