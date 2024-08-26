# This is equivalent to running the lua_inputs/config.yml file
# after building yodaStruct from seams-core
import bbdir.cyoda as cyoda

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

# Calculate the neighborlist by ID
nList = cyoda.neighListO(
    rcutoff=3.5,
    yCloud=resCloud,
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

# Hydrogen-bonded network using indices not IDs
hbnList = cyoda.neighbourListByIndex(
    yCloud=resCloud,
    nList=hbnList,
)

# Gets every ring (non-primitives included)
rings = cyoda.ringNetwork(
    nList=hbnList,
    maxDepth=6,
)

# Does the prism analysis for quasi-one-dimensional ice
cyoda.prismAnalysis(
    path="runOne/",  # outDir
    rings=rings,
    nList=hbnList,
    yCloud=resCloud,
    maxDepth=6,
    atomID=0,
    firstFrame=1,  # targetFrame
    currentFrame=1,  # frame
    doShapeMatching=False,
)
