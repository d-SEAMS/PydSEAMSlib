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
    pcd = _ase.to_pointcloud(atms)
    assert(pcd.box == resCloud.box)

    # Calculate the neighborlist by ID
    nList = cyoda.neighListO(
        rcutoff = 3.5,
        yCloud = resCloud,
        typeI = 2, #oxygenAtomType
    )

    verify(pprint.pformat(nList))
