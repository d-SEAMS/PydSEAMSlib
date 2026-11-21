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
