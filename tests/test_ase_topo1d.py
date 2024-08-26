import pprint
import shutil
import numpy as np

from approvaltests import verify, verify_file
from approvaltests.namer.default_namer_factory import NamerFactory
from ase.io import read as aseread


from pydseamslib import cyoda
from pydseamslib.adapters import _ase
from pathlib import Path

TRAJ = Path("subprojects/seams-core/input/traj/exampleTraj.lammpstrj")
strTRJ = str(TRAJ.absolute())


def test_nlist():
    # Get the frame
    resCloud = cyoda.readLammpsTrjreduced(
        filename=strTRJ,
        targetFrame=1,
        typeI=2,  # oxygenAtomType
        isSlice=False,
        coordLow=[0, 0, 0],
        coordHigh=[0, 0, 0],
    )

    # Construct a pointcloud
    atms = aseread(strTRJ)
    # TODO(ruhi): How should these be passed
    lammps_to_ase = {1: "H", 2: "O"}
    atms = _ase.map_LAMMPS_IDs_to_atomic_symbols(lammps_to_ase, atms)
    only_O_mask = [x.symbol == "O" for x in atms]
    molOID = np.repeat(np.arange(1, sum(only_O_mask) + 1), 1)
    pcd = _ase.to_pointcloud(
        atms, lammps_to_ase, only_O_mask, molOID, TRAJ, currentFrame=[1]
    )
    assert pcd.box == resCloud.box
    assert pcd.nop == resCloud.nop
    assert pcd.idIndexMap == resCloud.idIndexMap
    assert pcd.currentFrame == resCloud.currentFrame
    assert pcd.boxLow == resCloud.boxLow

    for idx in range(len(pcd.pts)):
        assert pcd.pts[idx].x == resCloud.pts[idx].x
        assert pcd.pts[idx].y == resCloud.pts[idx].y
        assert pcd.pts[idx].z == resCloud.pts[idx].z
        assert pcd.pts[idx].c_type == resCloud.pts[idx].c_type
        assert pcd.pts[idx].inSlice == resCloud.pts[idx].inSlice
        assert pcd.pts[idx].atomID == resCloud.pts[idx].atomID
        assert pcd.pts[idx].molID == resCloud.pts[idx].molID

    # Calculate the neighborlist by ID
    nl = cyoda.neighListO(
        rcutoff=3.5,
        yCloud=pcd,
        typeI=2,  # oxygenAtomType
    )

    # Validate the results
    verify(pprint.pformat(nl))


def test_hbnlist():
    # Get the frame, make a pointcloud
    atms = aseread(strTRJ)
    # TODO(ruhi): How should these be passed
    lammps_to_ase = {1: "H", 2: "O"}
    atms = _ase.map_LAMMPS_IDs_to_atomic_symbols(lammps_to_ase, atms)
    only_O_mask = [x.symbol == "O" for x in atms]
    molOID = np.repeat(np.arange(1, sum(only_O_mask) + 1), 1)
    pcd = _ase.to_pointcloud(
        atms, lammps_to_ase, only_O_mask, molOID, TRAJ, currentFrame=[1]
    )

    # Calculate the neighborlist by ID
    nl = cyoda.neighListO(
        rcutoff=3.5,
        yCloud=pcd,
        typeI=2,  # oxygenAtomType
    )

    # Get the hydrogen-bonded network for the current frame
    hl = cyoda.populateHbonds(
        filename=strTRJ,
        yCloud=pcd,
        nList=nl,
        targetFrame=1,
        Htype=1,  # hydrogen atom type
    )

    # Validate the results
    verify(pprint.pformat(hl))


def test_hbnlist1():
    # Get the frame, make a pointcloud
    atms = aseread(strTRJ)
    # TODO(ruhi): How should these be passed
    lammps_to_ase = {1: "H", 2: "O"}
    atms = _ase.map_LAMMPS_IDs_to_atomic_symbols(lammps_to_ase, atms)
    only_O_mask = [x.symbol == "O" for x in atms]
    molOID = np.repeat(np.arange(1, sum(only_O_mask) + 1), 1)
    pcd = _ase.to_pointcloud(
        atms, lammps_to_ase, only_O_mask, molOID, TRAJ, currentFrame=[1]
    )

    # Calculate the neighborlist by ID
    nl = cyoda.neighListO(
        rcutoff=3.5,
        yCloud=pcd,
        typeI=2,  # oxygenAtomType
    )

    # Get the hydrogen-bonded network for the current frame
    hl = cyoda.populateHbonds(
        filename=strTRJ,
        yCloud=pcd,
        nList=nl,
        targetFrame=1,
        Htype=1,  # hydrogen atom type
    )
    # Hydrogen-bonded network using indices not IDs
    hL = cyoda.neighbourListByIndex(
        yCloud=pcd,
        nList=hl,
    )

    # Validate the results
    verify(pprint.pformat(hL))


def test_rings():
    # Get the frame, make a pointcloud
    atms = aseread(strTRJ)
    # TODO(ruhi): How should these be passed
    lammps_to_ase = {1: "H", 2: "O"}
    atms = _ase.map_LAMMPS_IDs_to_atomic_symbols(lammps_to_ase, atms)
    only_O_mask = [x.symbol == "O" for x in atms]
    molOID = np.repeat(np.arange(1, sum(only_O_mask) + 1), 1)
    pcd = _ase.to_pointcloud(
        atms, lammps_to_ase, only_O_mask, molOID, TRAJ, currentFrame=[1]
    )

    # Calculate the neighborlist by ID
    nl = cyoda.neighListO(
        rcutoff=3.5,
        yCloud=pcd,
        typeI=2,  # oxygenAtomType
    )

    # Get the hydrogen-bonded network for the current frame
    hl = cyoda.populateHbonds(
        filename=strTRJ,
        yCloud=pcd,
        nList=nl,
        targetFrame=1,
        Htype=1,  # hydrogen atom type
    )
    # Hydrogen-bonded network using indices not IDs
    hL = cyoda.neighbourListByIndex(
        yCloud=pcd,
        nList=hl,
    )
    # Gets every ring (non-primitives included)
    Rgs = cyoda.ringNetwork(
        nList=hL,
        maxDepth=6,
    )

    # Validate the results
    verify(pprint.pformat(Rgs))


def test_prisms():
    # Construct a pointcloud
    atms = aseread(strTRJ)
    # TODO(ruhi): How should these be passed
    lammps_to_ase = {1: "H", 2: "O"}
    atms = _ase.map_LAMMPS_IDs_to_atomic_symbols(lammps_to_ase, atms)
    only_O_mask = [x.symbol == "O" for x in atms]
    molOID = np.repeat(np.arange(1, sum(only_O_mask) + 1), 1)
    pcd = _ase.to_pointcloud(
        atms, lammps_to_ase, only_O_mask, molOID, TRAJ, currentFrame=[1]
    )
    # Calculate the neighborlist by ID
    nl = cyoda.neighListO(
        rcutoff=3.5,
        yCloud=pcd,
        typeI=2,  # oxygenAtomType
    )
    # Get the hydrogen-bonded network for the current frame
    hl = cyoda.populateHbonds(
        filename=strTRJ,
        yCloud=pcd,
        nList=nl,
        targetFrame=1,
        Htype=1,  # hydrogen atom type
    )
    # Hydrogen-bonded network using indices not IDs
    hL = cyoda.neighbourListByIndex(
        yCloud=pcd,
        nList=hl,
    )
    # Gets every ring (non-primitives included)
    Rgs = cyoda.ringNetwork(
        nList=hL,
        maxDepth=6,
    )
    # Ensure the directory is not present before beginning
    outDir = "runOne/"  # / is important for the C++ engine..
    if Path(outDir).exists():
        shutil.rmtree(outDir)
    # Does the prism analysis for quasi-one-dimensional ice
    cyoda.prismAnalysis(
        path=outDir,
        rings=Rgs,
        nList=hL,
        yCloud=pcd,
        maxDepth=6,
        atomID=0,
        firstFrame=1,  # targetFrame
        currentFrame=1,  # frame
        doShapeMatching=False,
    )
    # Validate the run results
    verify_file(
        Path(f"{outDir}/topoINT/dataFiles/system-prisms-1.data"),
        options=NamerFactory.with_parameters("systemPrisms"),
    )
    verify_file(
        Path(f"{outDir}/topoINT/nPrisms.dat"),
        options=NamerFactory.with_parameters("nPrisms"),
    )
