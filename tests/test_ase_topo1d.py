import pprint

from approvaltests import verify

from pyseams import cyoda

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

    # Calculate the neighborlist by ID
    nList = cyoda.neighListO(
        rcutoff = 3.5,
        yCloud = resCloud,
        typeI = 2, #oxygenAtomType
    )
    verify(pprint.pformat(nList))
