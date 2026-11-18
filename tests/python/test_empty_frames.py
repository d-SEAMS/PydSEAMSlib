# Copyright (c) d-SEAMS developers
"""Empty frames on the pydseams surface.

Three inputs: no atoms, atoms that match no molecule or site kind, and
atoms that have no neighbour inside the cutoff. Each call returns an
empty array or list of the right shape, or a documented exception.
"""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms

import pydseams
import pydseamslib
from pydseams import band, features, yoda

BOX = (30.0, 30.0, 30.0)
ISOLATED = np.array(
    [[1.0, 1.0, 1.0], [12.0, 1.0, 1.0], [1.0, 12.0, 1.0]], dtype=np.float64
)
CUTOFF = 3.5
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


def empty_cloud(box=BOX):
    cloud = yoda.PointCloudDouble()
    cloud.nop = 0
    cloud.currentFrame = 1
    cloud.box = [float(v) for v in box]
    cloud.boxLow = [0.0, 0.0, 0.0]
    cloud.pts = []
    cloud.idIndexMap = {}
    return cloud


def isolated_frame():
    return pydseams.from_arrays(ISOLATED, BOX, cutoff=CUTOFF, bonded="cutoff")


def _rows(value):
    return [list(row) for row in value]


def test_constructors_reject_zero_atoms(tmp_path):
    with pytest.raises(ValueError, match="no atoms to load"):
        pydseams.from_arrays(np.zeros((0, 3)), BOX)
    with pytest.raises(ValueError, match="no atoms matched"):
        pydseams.from_ase(Atoms(cell=BOX, pbc=True), select=None)
    hydrogens = Atoms("H3", positions=ISOLATED, cell=BOX, pbc=True)
    with pytest.raises(ValueError, match="no atoms matched"):
        pydseams.from_ase(hydrogens, select="O")
    path = tmp_path / "empty.lammpstrj"
    path.write_text(EMPTY_DUMP, encoding="utf-8")
    with pytest.raises(ValueError, match="no atoms of type 1 or 2"):
        pydseams.read(path, bonded="cutoff")
    with pytest.raises(ValueError, match="no atoms to load"):
        pydseamslib.from_arrays(np.zeros((0, 3)), BOX)


def test_yoda_zero_atom_cloud_is_empty():
    cloud = empty_cloud()
    assert cloud.nop == 0
    assert _rows(yoda.neighListO(rcutoff=CUTOFF, yCloud=cloud, typeI=1)) == []
    assert _rows(yoda.neighList(rcutoff=CUTOFF, yCloud=cloud, typeI=1, typeJ=1)) == []
    assert _rows(yoda.halfNeighList(rcutoff=CUTOFF, yCloud=cloud, typeI=1)) == []
    assert _rows(yoda.getNewNeighbourListByIndex(cloud, CUTOFF)) == []
    assert _rows(yoda.kNearestNeighbourList(cloud, 4, CUTOFF + 1.5, 1, True)) == []
    mutual, union = yoda.kNearestNeighbourPair(cloud, 4, CUTOFF + 1.5, 1)
    assert _rows(mutual) == [] and _rows(union) == []
    assert _rows(yoda.neighbourListByIndex(yCloud=cloud, nList=[])) == []
    assert _rows(yoda.ringNetwork([], 6)) == []
    assert _rows(yoda.clearNeighbourList([[1, 2]])) == []
    assert list(yoda.mutualNearestUnlike(cloud, 1, 2)) == []
    assert _rows(yoda.voronoiFeatures(cloud, CUTOFF)) == []
    assert _rows(yoda.soapSpectrumAll(cloud, [], 3, 6, CUTOFF)) == []
    assert list(yoda.classifyTemplates(cloud, [], 12)) == []
    assert list(yoda.findBySignature([], [], "512")) == []
    hc, ddc = yoda.seededCageAffiliation([], [], [], [], False)
    assert list(hc) == [] and list(ddc) == []

    order = yoda.steinhardtQl(cloud, [], 6)
    assert list(order.ql) == [] and list(order.qlBar) == []
    voronoi = yoda.steinhardtQlVoronoi(cloud, CUTOFF, 6)
    assert list(voronoi.ql) == [] and list(voronoi.qlBar) == []
    yoda.getCorrelPlus(yCloud=cloud, nList=[], isSlice=False)
    yoda.getIceTypePlusNoPrint(yCloud=cloud, nList=[], isSlice=False)

    radius, values = yoda.partialRdf(yCloud=cloud, typeI=1, typeJ=1, rmax=6.0, nbins=12)
    assert len(radius) == 12 and len(values) == 12
    assert all(value == 0.0 for value in values)
    hist = yoda.partialRdfHist(cloud, 1, 1, 6.0, 12)
    assert hist.nI == 0 and hist.nJ == 0
    assert hist.volume == 30.0**3
    assert len(hist.r) == 12 and all(value == 0.0 for value in hist.g)

    density = yoda.densityZ(cloud, 0, 4, 2)
    assert list(density.z) == [3.75, 11.25, 18.75, 26.25]
    assert list(density.rho) == [0.0, 0.0, 0.0, 0.0]

    ions = yoda.ionEnvironment(cloud, [], [], 1, CUTOFF)
    assert list(ions.ion) == []
    assert ions.nIce == 0 and ions.nFront == 0 and ions.nLiquid == 0
    guests = yoda.guestOccupancy(cloud, [], [], 4.0)
    assert guests.occupied == 0 and guests.multiply == 0 and guests.free == 0
    assert list(guests.guestsPerCage) == [] and list(guests.cageOfGuest) == []
    fingerprint = yoda.topologyFingerprint([], 2, 7, [])
    assert list(fingerprint.atomKeys) == []
    assert list(fingerprint.ringCensus) == [0, 0, 0, 0, 0, 0, 0, 0]
    census = list(yoda.shellRingCensus([], [], 6))
    assert census == [0, 0, 0, 0, 0, 0, 0]
    domain = yoda.largestDomain(cloud, [], [])
    assert domain.largest == 0 and domain.subset == 0 and domain.percolation == 0.0
    typed = yoda.getPointCloudOneAtomType(
        cloud, empty_cloud(), 1, False, [0, 0, 0], [0, 0, 0]
    )
    assert typed.nop == 0
    assert yoda.moleculesInSingleSlice(cloud, True, [0, 0, 0], [1, 1, 1]) is None
    assert cloud.nop == 0
    edge = yoda.getEdgeMoleculesInRings(
        [], cloud, empty_cloud(), [0, 0, 0], [1, 1, 1], False
    )
    assert edge is None

    hydrogen = empty_cloud()
    bonded = yoda.populateHbondsWithInputClouds(yCloud=cloud, hCloud=hydrogen, nList=[])
    assert _rows(bonded) == []


def test_zero_molecules_leave_site_results_empty():
    frame = isolated_frame()
    table = yoda.parseSiteSpec("8=cationHead,9=anion")
    ions = frame.ion_cloud(table)
    assert ions.nop == 0
    pairs = frame.pairs(table)
    assert pairs.count == 0 and pairs.n_cation == 0 and pairs.n_anion == 0
    domain = frame.domain(table, "cationHead")
    assert domain.n == 0 and domain.largest == 0 and domain.percolation == 0.0
    with pytest.raises(ValueError, match="Hydrogen-bond"):
        list(frame.hbonds)
    donors = frame.hbonds_from_donors([], h_cloud=empty_cloud())
    assert _rows(donors) == [[1], [2], [3]]


def test_isolated_atoms_have_no_cutoff_neighbours(tmp_path):
    frame = isolated_frame()
    assert frame.n_atoms == 3
    assert _rows(frame.neighbor_list) == [[1], [2], [3]]
    assert _rows(frame.bonds_by_index) == [[0], [1], [2]]
    assert list(frame.rings) == []

    chill = frame.chill()
    plus = frame.chill_plus()
    assert chill.water == 3 and sum(chill.values()) == 3
    assert plus.water == 3 and sum(plus.values()) == 3

    cages = frame.cages(seeded=False)
    assert cages.n_ih == 0 and cages.n_ic == 0 and cages.n_water == 3
    assert cages.hc == [False, False, False]
    assert cages.ddc == [False, False, False]
    seeded = frame.cages(seeded=True)
    assert seeded.n_ih == 0 and seeded.n_ic == 0 and seeded.n_water == 3
    assert list(frame.cages_by_signature("512")) == []
    affiliation = frame.cage_affiliation()
    assert affiliation["six_rings"] == []
    assert affiliation["hc"] == [] and affiliation["ddc"] == []
    assert affiliation["reclassified"] == 0

    radius, values = frame.rdf(1, 1, cutoff=6.0, binwidth=0.5)
    assert list(radius) == [0.25 + 0.5 * i for i in range(12)]
    assert list(values) == [0.0] * 12
    assert frame.cn(1, 1, cutoff=3.0, binwidth=0.1) == 0.0
    running = frame.running_cn(1, 1, cutoff=3.0, binwidth=0.5)
    assert list(running) == [0.0] * 6

    profile = frame.density(bins=4, axis="z")
    dz = 30.0 / 4.0
    slab = (30.0 * 30.0) * dz
    assert profile.centres == tuple(dz * (i + 0.5) for i in range(4))
    assert profile.rho[0] == 3.0 / slab
    assert profile.rho[1:] == (0.0, 0.0, 0.0)

    order = frame.steinhardt(6)
    assert order["ql"] == [0.0, 0.0, 0.0]
    assert order["ql_bar"] == [0.0, 0.0, 0.0]
    weighted = frame.steinhardt_voronoi(6)
    assert len(weighted["ql"]) == 3 and len(weighted["ql_bar"]) == 3
    assert all(np.isfinite(weighted["ql"])) and all(np.isfinite(weighted["ql_bar"]))
    features_v = frame.voronoi_features()
    assert len(features_v) == 3 and all(
        len(row) == 3 and all(np.isfinite(row)) for row in features_v
    )

    spectra = frame.soap()
    assert len(spectra) == 3
    assert spectra[0] and all(value == 0.0 for row in spectra for value in row)
    hits = frame.classify_templates()
    assert len(hits) == 3
    assert all(hit["name"] == "other" and hit["kind"] == "other" for hit in hits)
    assert all(hit["rmsd"] == 1e300 for hit in hits)

    fingerprint = frame.fingerprint()
    assert len(fingerprint.atomKeys) == 3
    assert all(count == 0 for count in fingerprint.ringCensus)
    ions = frame.ion_environment([8])
    assert list(ions.ion) == [] and ions.nIce == 0
    shell = frame.hydration_shell_rings([8])
    assert shell[1] == [] or list(shell[1]) == []
    guests = frame.guest_occupancy([], [9])
    assert guests.occupied == 0 and guests.free == 0 and list(guests.cageOfGuest) == []
    both = frame.guest_occupancy_both("512", [9])
    assert both.radius.occupied == 0 and both.inside.occupied == 0

    with pytest.raises(ValueError, match="Hydrogen-bond"):
        frame.find_prisms(tmp_path / "prisms")
    monolayer = frame.monolayer_rings(tmp_path / "mono", 100.0, max_depth=4)
    assert monolayer[3]["count"] == 0 and monolayer[3]["coverage_xy"] == 0.0
    assert monolayer[4]["count"] == 0 and monolayer[4]["coverage_xy"] == 0.0
    centres, gram = frame.rdf_2d(tmp_path / "rdf2d", cutoff=6.0, binwidth=0.5)
    assert list(centres) == [0.25 + 0.5 * i for i in range(12)]
    assert len(gram) == 12 and all(np.isnan(value) for value in gram)


def test_featurizer_and_band_on_isolated_atoms():
    frame = isolated_frame()
    vector, states = features.IceFeaturizer(frame, chill=True).frame_features()
    assert vector.shape == (len(features.FEATURE_NAMES),)
    assert vector.dtype == np.float64
    assert states.shape == (3,) and states.dtype == np.int8
    assert np.all(states == features.STATE_WATER)
    water = features.FEATURE_NAMES.index("chill_water")
    assert vector[water] == 3.0
    assert np.all(vector[:water] == 0.0)
    assert np.all(vector[water + 1 :] == 0.0)
    rows = band.band_rows([frame])
    assert len(rows) == 1
    assert rows[0]["n_bonds"] == 0
    assert rows[0]["entered"] == [] and rows[0]["left"] == []
    assert rows[0]["rings"] == {}
    assert rows[0]["shortest_within_type"] == {1: 11.0}
