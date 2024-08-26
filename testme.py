import bbdir.cyoda as cyoda

trajectory = "subprojects/seams-core/input/traj/exampleTraj.lammpstrj"

# Get the frame
resCloud = cyoda.readLammpsTrjreduced(
    filename=trajectory,
    targetFrame=1,
    typeI=2,
    isSlice=False,
    coordLow=[0, 0, 0],
    coordHigh=[0, 0, 0],
)

# Calculate the neighborlist by ID
nList = cyoda.neighListO(3.5, resCloud, 2)

# Get the hydrogen-bonded network for the current frame
hbnList = cyoda.populateHbonds(trajectory, resCloud, nList, 1, 1)

# Hydrogen-bonded network using indices not IDs
hbnList = cyoda.neighbourListByIndex(yCloud=resCloud, nList=hbnList)

# Gets every ring (non-primitives included)
rings = cyoda.ringNetwork(hbnList, 6)

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
