"""Basic sanity checks for pydseams bindings."""

import pydseams
from pydseams import yoda


def test_version():
    assert hasattr(pydseams, "__version__")
    assert pydseams.__version__ == "2.10.0"


def test_import_yoda():
    assert hasattr(yoda, "PointDouble")
    assert hasattr(yoda, "PointCloudDouble")


def test_yoda_aliases():
    assert pydseams._core is pydseams.yoda
    assert pydseams.cyoda is pydseams.yoda


def test_point_double_construction():
    pt = yoda.PointDouble()
    pt.x = 1.0
    pt.y = 2.0
    pt.z = 3.0
    assert pt.x == 1.0
    assert pt.y == 2.0
    assert pt.z == 3.0


def test_pointcloud_double_construction():
    pcd = yoda.PointCloudDouble()
    pcd.nop = 0
    assert pcd.nop == 0


def test_readlammps_exists():
    assert hasattr(yoda, "readLammpsTrjreduced")


def test_neighlist_exists():
    assert hasattr(yoda, "neighListO")


def test_populate_hbonds_exists():
    assert hasattr(yoda, "populateHbonds")


def test_partial_rdf_exists():
    assert hasattr(yoda, "partialRdf")


def test_v240_symbols_exist():
    assert hasattr(yoda, "neighListPair")
    assert hasattr(yoda, "partialRdfHist")
    assert hasattr(yoda, "coordinationNumber")
    assert hasattr(yoda, "runningCN")
    assert hasattr(yoda, "firstMinimumBin")
    assert hasattr(yoda, "SiteTable")
    assert hasattr(yoda, "SiteKind")
    assert hasattr(yoda, "SiteFamily")
    assert yoda.Kind is yoda.SiteKind
    assert yoda.Family is yoda.SiteFamily
    assert hasattr(yoda, "parseSiteSpec")
    assert hasattr(yoda, "ionCloud")
    assert hasattr(yoda, "populateHbondsFromDonors")
    assert hasattr(yoda, "donatedHydrogenBond")


def test_ring_network_exists():
    assert hasattr(yoda, "ringNetwork")
    assert hasattr(yoda, "RingUpdater")


def test_lookup_table_q4_vec_length():
    result = yoda.lookupTableQ4Vec(angles=[0.7, 1.2])
    assert len(result) == 9


def test_lookup_table_q4_matches_vec():
    angles = [0.7, 1.2]
    vec = yoda.lookupTableQ4Vec(angles=angles)
    for m in range(9):
        single = yoda.lookupTableQ4(m=m, angles=angles)
        assert single == vec[m]
