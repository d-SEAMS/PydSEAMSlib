import pprint
import numpy as np

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
    molOID = np.repeat(np.arange(1,sum(only_O_mask)+1),1)
    pcd = _ase.to_pointcloud(atms,lammps_to_ase,only_O_mask,molOID)
    assert(pcd.box == resCloud.box)
    assert(pcd.nop == resCloud.nop)
    assert(pcd.idIndexMap == resCloud.idIndexMap)

    for idx in range(len(pcd.pts)):
        assert(pcd.pts[idx].x == resCloud.pts[idx].x)
        assert(pcd.pts[idx].y == resCloud.pts[idx].y)
        assert(pcd.pts[idx].z == resCloud.pts[idx].z)
        assert(pcd.pts[idx].c_type == resCloud.pts[idx].c_type)
        assert(pcd.pts[idx].inSlice == resCloud.pts[idx].inSlice)
        assert(pcd.pts[idx].atomID == resCloud.pts[idx].atomID)
        assert(pcd.pts[idx].molID == resCloud.pts[idx].molID)
    
   
    # Calculate the neighborlist by ID
    nList = cyoda.neighListO(
        rcutoff = 3.5,
        yCloud = resCloud,
        typeI = 2, #oxygenAtomType
    )
    nl = cyoda.neighListO(
        rcutoff = 3.5,
        yCloud = pcd,
        typeI = 2, #oxygenAtomType
    )
    for idx in range(len(nList)):
        assert(nList[idx]== nl[idx])

    verify(pprint.pformat(nList))


def test_hbnlist():
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
    molOID = np.repeat(np.arange(1,sum(only_O_mask)+1),1)
    pcd = _ase.to_pointcloud(atms,lammps_to_ase,only_O_mask,molOID)    

    # Calculate the neighborlist by ID
    nList = cyoda.neighListO(
        rcutoff = 3.5,
        yCloud = resCloud,
        typeI = 2, #oxygenAtomType
    )
    nl = cyoda.neighListO(
        rcutoff = 3.5,
        yCloud = pcd,
        typeI = 2, #oxygenAtomType
    )

    #Get the hydrogen-bonded network for the current frame
    hbnList = cyoda.populateHbonds(
        filename = trajectory,
        yCloud = resCloud,
        nList = nList, 
        targetFrame = 1,
        Htype = 1, #hydrogen atom type
    )
    hl = cyoda.populateHbonds(
        filename = trajectory,
        yCloud = pcd,
        nList = nl, 
        targetFrame = 1,
        Htype = 1, #hydrogen atom type
    )
    for idx in range(len(hbnList)):
        assert(hbnList[idx]== hl[idx])


    verify(pprint.pformat(hl))

def test_hbnlist1():
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
    molOID = np.repeat(np.arange(1,sum(only_O_mask)+1),1)
    pcd = _ase.to_pointcloud(atms,lammps_to_ase,only_O_mask,molOID)    

    # Calculate the neighborlist by ID
    nList = cyoda.neighListO(
        rcutoff = 3.5,
        yCloud = resCloud,
        typeI = 2, #oxygenAtomType
    )
    nl = cyoda.neighListO(
        rcutoff = 3.5,
        yCloud = pcd,
        typeI = 2, #oxygenAtomType
    )

    #Get the hydrogen-bonded network for the current frame
    hbnList = cyoda.populateHbonds(
        filename = trajectory,
        yCloud = resCloud,
        nList = nList, 
        targetFrame = 1,
        Htype = 1, #hydrogen atom type
    )
    hl = cyoda.populateHbonds(
        filename = trajectory,
        yCloud = pcd,
        nList = nl, 
        targetFrame = 1,
        Htype = 1, #hydrogen atom type
    )
    #Hydrogen-bonded network using indices not IDs
    hbnList =  cyoda.neighbourListByIndex(
        yCloud = resCloud,
        nList = hbnList,
    )
    hL =  cyoda.neighbourListByIndex(
        yCloud = pcd,
        nList = hl,
    )
    for idx in range(len(hbnList)):
        assert(hbnList[idx]== hL[idx])


    verify(pprint.pformat(hbnList))
