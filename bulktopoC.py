import bbdir.cyoda as cyoda

trajectory = "subprojects/seams-core/input/traj/cluster-417.lammpstrj"

# Get the frame
resCloud = cyoda.readLammpsTrjreduced(
    filename=trajectory,
    targetFrame=1302,
    typeI=1,  # oxygenAtomType
    isSlice=False,
    coordLow=[0, 0, 0],
    coordHigh=[0, 0, 0],
)

# Calculate the neighborlist by ID
nList = cyoda.neighListO(
    rcutoff=3.5,
    yCloud=resCloud,
    typeI=1,  # oxygenAtomType
)
solCloud = cyoda.PointCloudDouble()
iceList = []
clump = cyoda.clusterAnalysis(
    path="runOne/",  # outDir
    iceCloud=solCloud,
    yCloud=resCloud,
    nList=nList,
    iceNeighbourList=iceList,
    cutoff=3.5,
    firstFrame=1302,
    bopAnalysis="q6",
)

# Gets every ring (non-primitives included)
rings = cyoda.ringNetwork(
    nList=iceList,
    maxDepth=6,
)

# Finds DDCs and HCs
tum3 = cyoda.topoUnitMatchingBulk(
    path="runOne/",  # outDir
    rings=rings,
    iceNeighbourList=iceList,
    yCloud=solCloud,
    firstFrame=1302,
    printClusters=True,
    onlyTetrahedral=False,
)
