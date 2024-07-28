import pprint

from approvaltests import verify
from ase.io import read as aseread


from pyseams import cyoda
from pyseams.adapters import _ase

def test_nlist():
    trajectory ="subprojects/seams-core/input/traj/exampleTraj.lammpstrj"
    # Get the frame
    resCloud = cyoda.readLammpsTrjreduced(
              filename = trajectory,
              targetFrame = 1,
              typeI = 2, #oxygenAtomType
              isSlice = False,
              coordLow = [0,0,0],
              coordHigh = [0,0,0],
    )

    # Construct a pointcloud
    atms = aseread(trajectory)
    # TODO(ruhi): How should these be passed
    lammps_to_ase = {1: 'H', 2: 'O'}
    atms = _ase.map_2(lammps_to_ase, atms)
    only_O_mask = [x.symbol == 'O' for x in atms]
    pcd = _ase.to_pointcloud(atms[only_O_mask],lammps_to_ase)
    assert(pcd.box == resCloud.box)
    assert(pcd.nop == resCloud.nop)
    assert(pcd.pts[0].x == resCloud.pts[0].x)
    assert(pcd.pts[0].y == resCloud.pts[0].y)
    assert(pcd.pts[0].z == resCloud.pts[0].z)
    assert(pcd.pts[0].c_type == resCloud.pts[0].c_type)


    # Calculate the neighborlist by ID
    nList = cyoda.neighListO(
        rcutoff = 3.5,
        yCloud = resCloud,
        typeI = 2, #oxygenAtomType
    )

    verify(pprint.pformat(nList))
