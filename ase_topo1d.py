from pydseamslib import cyoda
from pydseamslib.adapters import _ase
import numpy as np

from ase.io import read as aseread

trajectory = "subprojects/seams-core/input/traj/exampleTraj.lammpstrj"


# Get the frame
resCloud = cyoda.readLammpsTrjreduced(
    filename=trajectory,
    targetFrame=1,
    typeI=2,  # oxygenAtomType
    isSlice=False,
    coordLow=[0, 0, 0],
    coordHigh=[0, 0, 0],
)
# pprint.pprint(dir(resCloud))

atms = aseread(trajectory)
# In ASE, we want to work with atomic symbols instead of LAMMPS types
lammps_to_ase = {1: "H", 2: "O"}
atms = _ase.map_LAMMPS_IDs_to_atomic_symbols(lammps_to_ase, atms)
openfile = "subprojects/seams-core/input/traj/exampleTraj.lammpstrj"
only_O_mask = [x.symbol == "O" for x in atms]
# user has to provide proper molID for inslice, each molecule must have one molID. for eg: In H2O, both H atoms and O atom must have same molID.
# if one wants molHID, they can just change the last ,1) as ,2) in  the following expression:
molOID = np.repeat(np.arange(1, sum(only_O_mask) + 1), 1)
pcd = _ase.to_pointcloud(
    atms, lammps_to_ase, only_O_mask, molOID, openfile, currentFrame=[1]
)

# Calculate the neighborlist by ID
nList = cyoda.neighListO(
    rcutoff=3.5,
    yCloud=resCloud,
    typeI=2,  # oxygenAtomType
)
nl = cyoda.neighListO(
    rcutoff=3.5,
    yCloud=pcd,
    typeI=2,  # oxygenAtomType
)

# Get the hydrogen-bonded network for the current frame
hbnList = cyoda.populateHbonds(
    filename=trajectory,
    yCloud=resCloud,
    nList=nList,
    targetFrame=1,
    Htype=1,  # hydrogen atom type
)
hl = cyoda.populateHbonds(
    filename=trajectory,
    yCloud=pcd,
    nList=nl,
    targetFrame=1,
    Htype=1,  # hydrogen atom type
)
# Hydrogen-bonded network using indices not IDs
hbnList = cyoda.neighbourListByIndex(
    yCloud=resCloud,
    nList=hbnList,
)
hL = cyoda.neighbourListByIndex(
    yCloud=pcd,
    nList=hl,
)
# Gets every ring (non-primitives included)
rings = cyoda.ringNetwork(
    nList=hbnList,
    maxDepth=6,
)
Rgs = cyoda.ringNetwork(
    nList=hL,
    maxDepth=6,
)
# Does the prism analysis for quasi-one-dimensional ice
cyoda.prismAnalysis(
    path="runOne/",  # outDir
    rings=Rgs,
    nList=hL,
    yCloud=resCloud,
    maxDepth=6,
    atomID=0,
    firstFrame=1,  # targetFrame
    currentFrame=1,  # frame
    doShapeMatching=False,
)

# print(pprint.pformat(nList))
