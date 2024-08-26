import bbdir.cyoda as cyoda

trajectory = "subprojects/seams-core/input/traj/clathrate-thf.lammpstrj"

# Get the frame
resCloud = cyoda.readLammpsTrjreduced(
    filename=trajectory,
    targetFrame=1,
    typeI=1,  # oxygenAtomType
    isSlice=True,
    coordLow=[0, 0, 0],
    coordHigh=[34.728, 0, 0],
)

# Calculate the neighborlist by ID
nList = cyoda.neighListO(
    rcutoff=3.5,
    yCloud=resCloud,
    typeI=1,  # oxygenAtomType
)

# Get the hydrogen-bonded network for the current frame
hbnList = cyoda.populateHbonds(
    filename=trajectory,
    yCloud=resCloud,
    nList=nList,
    targetFrame=1,
    Htype=2,  # hydrogen atom type
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

# Writes out primitive rings for a bulk system
ring = cyoda.bulkPolygonRingAnalysis(
    path="runOne/",  # outDir
    rings=rings,
    nList=hbnList,
    yCloud=resCloud,
    maxDepth=6,
    firstFrame=1,
)
