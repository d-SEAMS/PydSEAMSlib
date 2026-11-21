import pytest

from pydseams import iter_trajectory, read


def test_topologized_water_sites(tmp_path):
    sites = [(1, "WAT", "O", .50), (1, "WAT", "H1", .59),
             (2, "WAT", "O", .78), (2, "WAT", "H1", .69),
             (3, "ION", "NA", .90)]
    gro = tmp_path / "sites.gro"
    with gro.open("w") as f:
        for step in range(2):
            f.write(f"water t= {step}.0\n{len(sites)}\n")
            for i, (resid, resname, name, x) in enumerate(sites, 1):
                f.write(f"{resid:5d}{resname:<5}{name:>5}{i:5d}{x+step*.01:8.3f}{.5:8.3f}{.5:8.3f}\n")
            f.write("   2.00000   2.00000   2.00000\n")
    options = dict(atom_name="O", residue_name="WAT", hydrogen_names=("H1",))
    frames = list(iter_trajectory(gro, **options))
    assert len(frames) == 2
    assert frames[0].positions[0][0] == pytest.approx(5.0)
    assert frames[1].positions[0][0] == pytest.approx(5.1)
    assert frames[0].n_atoms == 2
    assert [p.molID for p in frames[0].cloud.pts] == [1, 2]
    assert 2 in list(frames[0].hbonds[0])[1:]
    assert read(gro, frame=2, **options).positions == frames[1].positions


def test_xtc_selection_preserves_tilt_and_frame_numbers(tmp_path):
    import struct
    import numpy as np

    # XTC stores at most nine sites as uncompressed XDR coordinate triples.
    # Both off-diagonal terms and the raw atom ordering must survive selection.
    cell = [2., 0., 0., .4, 2., 0., .1, .2, 2.]
    sites = [("SOLX", "OW"), ("SOL", "OW"), ("SOL", "HW1"), ("SOL", "HW2")]
    gro = tmp_path / "water.gro"
    records = [f"{1 if r == 'SOLX' else 2:5d}{r:<5}{a:>5}{i:5d}{.2:8.3f}{.2:8.3f}{.2:8.3f}"
               for i, (r, a) in enumerate(sites, 1)]
    gro.write_text("\n".join(["labelled sites", "4", *records, "2 2 2 0 0 .4 0 .1 .2", ""]))
    xtc = tmp_path / "water.xtc"
    with xtc.open("wb") as stream:
        for i in range(68):
            xyz = [.9, .9, .9, .5+i*.001, .5, .5, .6, .5, .5, .5, .6, .5]
            stream.write(struct.pack(">3if9fi12f", 1995, 4, i, i*.5, *cell, 4, *xyz))
    opts = dict(topology=gro, residue_name="SOL", atom_name="OW", hydrogen_names=("HW1", "HW2"))
    frames = list(iter_trajectory(xtc, start=1, stop=68, stride=2, **opts))
    assert [f.frame for f in frames] == list(range(2, 69, 2))
    assert all(f.n_atoms == 1 and f._h_cloud.nop == 2 for f in frames)
    assert np.allclose(frames[-1].positions, [[5.67, 5., 5.]])
    assert np.allclose(list(frames[-1].cloud.box)[3:], [4., 1., 2.])


def test_incomplete_named_hydrogens(tmp_path):
    sites = [(1, "O"), (1, "H1"), (1, "H2"), (2, "O"), (2, "H1")]
    path = tmp_path / "incomplete.gro"
    atoms = [f"{resid:5d}{'WAT':<5}{name:>5}{i:5d}{.5:8.3f}{.5:8.3f}{.5:8.3f}"
             for i, (resid, name) in enumerate(sites, 1)]
    path.write_text("\n".join(["incomplete WAT", "5", *atoms, "2 2 2", ""]))
    with pytest.raises(ValueError, match="H2"):
        list(iter_trajectory(path, atom_name="O", hydrogen_names=("H1", "H2")))
