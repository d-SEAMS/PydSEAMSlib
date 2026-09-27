"""Public Frame API: no scratch files, cages, ASE optional."""

from pathlib import Path

import pytest

from pydseams import (
    CageScore,
    ContactPairs,
    DensityProfile,
    DomainStats,
    Frame,
    IceCounts,
    Trajectory,
    from_arrays,
    read,
)

TRAJ = Path(__file__).resolve().parents[1] / "data" / "exampleTraj.lammpstrj"


def test_read_guesses_oxygen():
    frame = read(TRAJ)
    assert frame.n_atoms == 250
    assert frame.atom_type == 2


def test_read_all_atoms_keeps_mixed_lammps_types():
    frame = read(TRAJ, all_atoms=True, atom_type=2)
    assert frame.n_atoms == 750
    assert frame.atom_type == 2
    assert {point.c_type for point in frame.cloud.pts} == {1, 2}


def test_read_all_atoms_rejects_coordinate_slice():
    with pytest.raises(ValueError, match="all_atoms and region"):
        read(
            TRAJ,
            all_atoms=True,
            region=([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        )


def test_chill_plus_no_tempdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    frame = read(TRAJ.resolve())
    counts = frame.chill_plus()
    assert isinstance(counts, IceCounts)
    assert counts.interClathrate == 12
    assert counts.water == 238
    assert counts.cubic == 0
    assert counts.hexagonal == 0
    assert sum(counts.values()) == 250
    leftover = [p for p in tmp_path.iterdir() if p.name.startswith("dseams_")]
    assert leftover == []


def test_chill_plus_replaces_chill_classification():
    frame = read(TRAJ)
    frame.chill()
    counts = frame.chill_plus()
    assert counts.interClathrate == 12
    assert counts.water == 238
    assert sum(counts.values()) == frame.n_atoms


def test_find_by_signature_hex_prism():
    from pydseams import yoda

    rings = [
        [0, 1, 2, 3, 4, 5],
        [6, 7, 8, 9, 10, 11],
        [0, 1, 7, 6],
        [1, 2, 8, 7],
        [2, 3, 9, 8],
        [3, 4, 10, 9],
        [4, 5, 11, 10],
        [5, 0, 6, 11],
    ]
    n_list = [[] for _ in range(12)]
    found = yoda.findBySignature(rings, n_list, "4:6,6:2")
    assert len(found) == 1
    assert found[0]["signature"] == "4:6,6:2"
    assert len(found[0]["vertices"]) == 12
    assert len(found[0]["faces"]) == 8


def test_cages_by_signature_on_mixed_water():
    frame = read(TRAJ, bonded="cutoff")
    found = frame.cages_by_signature("hc")
    assert isinstance(found, list)


def test_cages_seeded_on_mixed_water():
    frame = read(TRAJ, bonded="cutoff")
    score = frame.cages(seeded=True)
    assert isinstance(score, CageScore)
    assert score.n_ih == 0
    assert score.n_ic == 0
    assert score.n_water == frame.n_atoms


def test_trajectory_alias():
    assert Trajectory is Frame
    traj = Trajectory(str(TRAJ.resolve()))
    counts = traj.classify_chill_plus()
    assert counts["hexagonal"] == 0
    assert counts["water"] == 238
    assert counts["interClathrate"] == 12


def test_from_arrays_roundtrip_positions():
    src = read(TRAJ)
    frame = from_arrays(src.positions, src.box, numbers=[8] * src.n_atoms)
    assert frame.n_atoms == 250
    assert frame.positions[0] == src.positions[0]


def test_frame_rdf_two_type(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    positions = [
        [0.0, 0.0, 0.0],
        [2.8, 0.0, 0.0],
        [0.0, 2.8, 0.0],
        [2.8, 2.8, 0.0],
    ]
    frame = from_arrays(positions, [30.0, 30.0, 30.0], numbers=[1, 2, 1, 2])
    r, g = frame.rdf(1, 2)
    nbin = int(12.0 / 0.05)
    assert len(r) == nbin
    assert len(g) == nbin
    assert any(gi > 0 for gi in g)
    leftover = [p for p in tmp_path.iterdir() if p.name.startswith("dseams_")]
    assert leftover == []


def test_frame_cn_matches_unlike_degree():
    positions = [
        [0.0, 0.0, 0.0],
        [2.8, 0.0, 0.0],
        [0.0, 2.8, 0.0],
        [2.8, 2.8, 0.0],
    ]
    frame = from_arrays(positions, [30.0, 30.0, 30.0], numbers=[1, 2, 1, 2])
    assert frame.cn(1, 2, cutoff=3.0, binwidth=0.1) == pytest.approx(1.0)
    assert frame.cn(1, 2, cutoff=1.0, binwidth=0.1) == pytest.approx(0.0)


def test_frame_running_cn_last_bin_matches_cn():
    from pydseams import yoda

    positions = [
        [0.0, 0.0, 0.0],
        [2.8, 0.0, 0.0],
        [0.0, 2.8, 0.0],
        [2.8, 2.8, 0.0],
    ]
    frame = from_arrays(positions, [30.0, 30.0, 30.0], numbers=[1, 2, 1, 2])
    cn_run = frame.running_cn(1, 2, cutoff=3.0, binwidth=0.1)
    assert isinstance(cn_run, list)
    assert len(cn_run) == 30
    assert cn_run[-1] == pytest.approx(1.0)
    assert cn_run[0] == pytest.approx(0.0)
    assert frame.running_cn(1, 2, cutoff=1.0, binwidth=0.1)[-1] == pytest.approx(0.0)
    hist = yoda.partialRdfHist(yCloud=frame.cloud, typeI=1, typeJ=2, rmax=3.0, nbins=30)
    assert hist.nJ == 2
    assert hist.volume == pytest.approx(30.0**3)
    assert yoda.runningCN(h=hist, rhoJ=hist.nJ / hist.volume)[-1] == pytest.approx(1.0)


def test_site_table_type_one_is_not_chemistry():
    from pydseams import yoda

    table = yoda.parseSiteSpec("2=anion")
    assert table.ofType(1) == yoda.SiteKind.unspecified
    assert table.ofType(2) == yoda.SiteKind.anion


def test_ion_cloud_two_atom_com():
    from pydseams import yoda

    positions = [[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]]
    frame = from_arrays(positions, [10.0, 10.0, 10.0], numbers=[1, 1])
    frame.cloud.pts[0].molID = 7
    frame.cloud.pts[0].atomID = 1
    frame.cloud.pts[1].molID = 7
    frame.cloud.pts[1].atomID = 2
    table = yoda.parseSiteSpec("1=cationHead")
    ions = frame.ion_cloud(table)
    assert ions.nop == 1
    assert ions.pts[0].c_type == 1
    assert ions.pts[0].x == pytest.approx(0.0)


def test_ion_cloud_two_type_ints():
    from pydseams import yoda

    positions = [[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]]
    frame = from_arrays(positions, [10.0, 10.0, 10.0], numbers=[1, 2])
    frame.cloud.pts[0].molID = 7
    frame.cloud.pts[0].atomID = 1
    frame.cloud.pts[1].molID = 8
    frame.cloud.pts[1].atomID = 2
    ions = yoda.ionCloud(src=frame.cloud, cationType=1, anionType=2)
    assert ions.nop == 2
    assert sorted(p.c_type for p in ions.pts) == [1, 2]


def test_ion_cloud_type_to_kind_dict():
    from pydseams import yoda

    positions = [[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]]
    frame = from_arrays(positions, [10.0, 10.0, 10.0], numbers=[1, 2])
    frame.cloud.pts[0].molID = 7
    frame.cloud.pts[0].atomID = 1
    frame.cloud.pts[1].molID = 8
    frame.cloud.pts[1].atomID = 2
    ions = yoda.ionCloud(
        src=frame.cloud,
        typeToKind={1: yoda.Kind.cationHead, 2: yoda.Kind.anion},
    )
    assert ions.nop == 2
    assert sorted(p.c_type for p in ions.pts) == [1, 2]


def test_density_profile_by_type():
    positions = [
        [1.0, 1.0, 1.0],
        [2.0, 2.0, 3.0],
        [3.0, 3.0, 7.0],
    ]
    frame = from_arrays(positions, [10.0, 10.0, 10.0], numbers=[1, 1, 2])
    profile = frame.density(bins=2, axis="z", atom_type=1)
    assert isinstance(profile, DensityProfile)
    assert profile.axis == "z"
    assert profile.atom_type == 1
    assert profile.centres == pytest.approx((2.5, 7.5))
    assert profile.rho == pytest.approx((0.004, 0.0))


def test_density_profile_by_site_kind():
    from pydseams import yoda

    positions = [
        [1.0, 1.0, 1.0],
        [2.0, 2.0, 3.0],
        [3.0, 3.0, 7.0],
    ]
    frame = from_arrays(positions, [10.0, 10.0, 10.0], numbers=[1, 1, 2])
    table = yoda.parseSiteSpec("1=polar,2=apolar")
    profile = frame.density(
        bins=2,
        axis="z",
        table=table,
        kind=yoda.Kind.polar,
    )
    assert isinstance(profile, DensityProfile)
    assert profile.site_kind == "polar"
    assert profile.centres == pytest.approx((2.5, 7.5))
    assert profile.rho == pytest.approx((0.004, 0.0))


def test_contact_pairs_match_mutual_nearest_unlike():
    from pydseams import yoda

    positions = [
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [7.0, 0.0, 0.0],
        [8.0, 0.0, 0.0],
    ]
    frame = from_arrays(positions, [10.0, 10.0, 10.0], numbers=[1, 2, 2, 1])
    table = yoda.parseSiteSpec("1=cationHead,2=anion")
    result = frame.pairs(table)
    assert isinstance(result, ContactPairs)
    assert result.count == 2
    assert result.n_cation == 2
    assert result.n_anion == 2
    assert result.pairs == ((0, 1), (3, 2))


def test_largest_site_domain():
    from pydseams import yoda

    positions = [
        [1.0, 1.0, 1.0],
        [2.0, 1.0, 1.0],
        [3.0, 1.0, 1.0],
        [12.0, 1.0, 1.0],
        [15.0, 1.0, 1.0],
    ]
    frame = from_arrays(positions, [20.0, 20.0, 20.0], numbers=[1, 1, 1, 1, 2])
    table = yoda.parseSiteSpec("1=polar,2=apolar")
    result = frame.domain(table, yoda.Kind.polar, cutoff=1.1)
    assert isinstance(result, DomainStats)
    assert result.site_kind == "polar"
    assert result.n == 4
    assert result.largest == 3
    assert result.percolation == pytest.approx(0.75)


def test_low_level_workflow_bindings():
    from pydseams import yoda

    frame = from_arrays(
        [[1.0, 0.0, 1.0], [2.0, 0.0, 2.0]],
        [10.0, 10.0, 10.0],
        numbers=[1, 2],
    )
    density = yoda.densityZ(frame.cloud, 0, 2, 2)
    assert isinstance(density, yoda.DensityZ)
    assert sum(density.rho) == pytest.approx(0.004)
    assert yoda.mutualNearestUnlike(frame.cloud, 1, 2) == [(0, 1)]

    by_index = yoda.getNewNeighbourListByIndex(frame.cloud, 2.0)
    by_id = [[frame.cloud.pts[index].atomID for index in row] for row in by_index]
    domain = yoda.largestDomain(frame.cloud, by_id, [True, True])
    assert isinstance(domain, yoda.Domain)
    assert domain.subset == 2
    assert domain.largest == 2
    assert domain.percolation == pytest.approx(1.0)


def test_available_readers():
    from pydseams import available_readers

    readers = available_readers()
    assert readers["lammps"] is True
    assert "xyz" in readers


def test_read_xyz(tmp_path):
    xyz = tmp_path / "pair.xyz"
    xyz.write_text("2\n\nO 0.0 0.0 0.0\nO 1.0 0.0 0.0\n")
    frame = read(xyz)
    assert frame.n_atoms == 2


def test_from_xyz_keeps_a_nitrogen_pair(tmp_path):
    from pydseams import Frame

    xyz = tmp_path / "sin.xyz"
    xyz.write_text(
        "4\n\n"
        "Si 0.0 0.0 0.0\n"
        "Si 3.0 0.0 0.0\n"
        "N 0.0 3.0 0.0\n"
        "N 1.1 3.0 0.0\n"
    )
    frame = Frame.from_xyz(xyz, cutoff=2.0)
    types = {int(pt.c_type) for pt in frame.cloud.pts}
    assert types == {1, 2}
    # N is type 2. Its row must contain the other N, not only silicons.
    rows = frame.bonds_by_index
    nitrogen = [i for i, pt in enumerate(frame.cloud.pts) if int(pt.c_type) == 2]
    assert len(nitrogen) == 2
    assert nitrogen[1] in rows[nitrogen[0]]
    assert nitrogen[0] in rows[nitrogen[1]]


def test_to_solvis_optional():
    pytest.importorskip("solvis")
    from pydseams import to_solvis

    frame = read(TRAJ, bonded="cutoff")
    system = to_solvis(frame)
    assert len(system.atoms) == frame.n_atoms


def test_from_ase_optional():
    np = pytest.importorskip("numpy")
    pytest.importorskip("ase")
    from ase import Atoms
    from pydseams import from_ase

    src = read(TRAJ)
    atoms = Atoms(
        symbols=["O"] * src.n_atoms,
        positions=src.positions,
        cell=src.box,
        pbc=True,
    )
    frame = from_ase(atoms, select="O", bonded="cutoff")
    assert frame.n_atoms == 250
    back = frame.to_ase()
    assert len(back) == 250
    np.testing.assert_allclose(back.cell.array, atoms.cell.array, atol=1.0e-12)
    np.testing.assert_allclose(back.positions, atoms.positions, atol=1.0e-12)
    assert back.get_chemical_symbols() == atoms.get_chemical_symbols()


def test_from_ase_auto_assigns_hydrogens_to_nearest_oxygen():
    pytest.importorskip("ase")
    from ase import Atoms
    from pydseams import from_ase

    atoms = Atoms(
        "OHHOHH",
        positions=[
            [5.0, 5.0, 5.0],
            [4.0, 5.0, 5.0],
            [5.0, 4.0, 5.0],
            [7.8, 5.0, 5.0],
            [6.8, 5.0, 5.0],
            [7.8, 6.0, 5.0],
        ],
        cell=[20.0, 20.0, 20.0],
        pbc=True,
    )

    frame = from_ase(atoms)

    assert frame.bonded == "hbond"
    assert frame.hbonds == [[1, 2], [2, 1]]


def test_from_ase_preserves_mol_id_for_hydrogen_donors():
    np = pytest.importorskip("numpy")
    pytest.importorskip("ase")
    from ase import Atoms
    from pydseams import from_ase

    atoms = Atoms(
        "OHHOHH",
        positions=[
            [5.0, 5.0, 5.0],
            [4.0, 5.0, 5.0],
            [5.0, 4.0, 5.0],
            [7.8, 5.0, 5.0],
            [6.8, 5.0, 5.0],
            [7.8, 6.0, 5.0],
        ],
        cell=[20.0, 20.0, 20.0],
        pbc=True,
    )
    atoms.new_array("mol-id", np.array([10, 10, 10, 20, 20, 20]))

    frame = from_ase(atoms)

    assert [point.molID for point in frame.cloud.pts] == [10, 20]
    assert [point.molID for point in frame._h_cloud.pts] == [10, 10, 20, 20]


def test_from_ase_all_atom_auto_uses_cutoff_bonding():
    pytest.importorskip("ase")
    from ase import Atoms
    from pydseams import from_ase

    atoms = Atoms(
        "OH",
        positions=[[1.0, 1.0, 1.0], [1.9, 1.0, 1.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=True,
    )

    frame = from_ase(atoms, select=None)

    assert frame.bonded == "cutoff"


def test_from_ase_hbond_rejects_selected_hydrogens():
    pytest.importorskip("ase")
    from ase import Atoms
    from pydseams import from_ase

    atoms = Atoms(
        "OH",
        positions=[[1.0, 1.0, 1.0], [1.9, 1.0, 1.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=True,
    )

    with pytest.raises(ValueError, match="heavy-atom selection"):
        from_ase(atoms, select=None, bonded="hbond")


def test_to_ase_exports_chill_plus_labels():
    pytest.importorskip("ase")
    from ase import Atoms
    from pydseams import from_ase

    atoms = Atoms(
        "OOOO",
        positions=[
            [1.0, 1.0, 1.0],
            [2.8, 1.0, 1.0],
            [1.0, 2.8, 1.0],
            [2.8, 2.8, 1.0],
        ],
        cell=[12.0, 12.0, 12.0],
        pbc=True,
    )
    frame = from_ase(atoms, bonded="cutoff")

    frame.chill_plus()
    labelled = frame.to_ase()

    assert list(labelled.arrays["ice_type"]) == [
        point.iceType.name for point in frame.cloud.pts
    ]


def test_to_ase_exports_cage_flags():
    pytest.importorskip("ase")
    from ase import Atoms
    from pydseams import from_ase

    atoms = Atoms(
        "OOOO",
        positions=[
            [1.0, 1.0, 1.0],
            [2.8, 1.0, 1.0],
            [1.0, 2.8, 1.0],
            [2.8, 2.8, 1.0],
        ],
        cell=[12.0, 12.0, 12.0],
        pbc=True,
    )
    frame = from_ase(atoms, bonded="cutoff")

    score = frame.cages(seeded=False)
    labelled = frame.to_ase()

    assert list(labelled.arrays["hc"]) == score.hc
    assert list(labelled.arrays["ddc"]) == score.ddc


def test_from_ase_roundtrips_general_periodic_cell():
    np = pytest.importorskip("numpy")
    pytest.importorskip("ase")
    from ase import Atoms
    from pydseams import from_ase

    angle = np.deg2rad(31.0)
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    restricted = np.array(
        [
            [4.0, 0.0, 0.0],
            [1.0, 5.0, 0.0],
            [0.5, 1.5, 6.0],
        ]
    )
    cell = restricted @ rotation
    scaled = np.array([[0.1, 0.2, 0.3], [0.7, 0.6, 0.5]])
    atoms = Atoms(
        symbols=["O", "Na"],
        scaled_positions=scaled,
        cell=cell,
        pbc=True,
    )
    atoms.set_celldisp([0.4, -0.3, 0.2])

    frame = from_ase(atoms, select=None, bonded="cutoff")
    assert len(frame.box) == 6
    back = frame.to_ase()
    np.testing.assert_allclose(back.cell.array, cell, atol=1.0e-12)
    np.testing.assert_allclose(back.positions, atoms.positions, atol=1.0e-12)
    np.testing.assert_allclose(back.get_celldisp(), atoms.get_celldisp(), atol=1.0e-12)
    assert back.get_chemical_symbols() == atoms.get_chemical_symbols()


@pytest.mark.parametrize(
    "cell,pbc,match",
    [
        ([10.0, 10.0, 10.0], [True, True, False], "periodic in all three"),
        ([[10.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 10.0]], True, "singular"),
    ],
)
def test_from_ase_rejects_unsupported_periodicity(cell, pbc, match):
    pytest.importorskip("ase")
    from ase import Atoms
    from pydseams import from_ase

    atoms = Atoms("O", positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=pbc)
    with pytest.raises(ValueError, match=match):
        from_ase(atoms, bonded="cutoff")
