import bbdir.cyoda as cyoda

trajectory = "subprojects/seams-core/input/traj/mW_cubic.lammpstrj"

solCloud = cyoda.PointCloudDouble()
iceList = []

# Get the frame
resCloud = cyoda.readLammpsTrjO(
    filename=trajectory,
    targetFrame=1,
    typeO=1,  # oxygenAtomType
    isSlice=False,
    coordLow=[0, 0, 0],
    coordHigh=[50, 0, 0],
)

# Calculate the neighborlist by ID
nList = cyoda.neighListO(
    rcutoff=3.5,
    yCloud=resCloud,
    typeI=1,  # oxygenAtomType
)

# Calculate Cij (cloud,slice)
resCloud = cyoda.getCorrelPlus(
    yCloud=resCloud,
    nList=nList,
    isSlice=False,
)

# Write out data (cloud,slice,name)
resCloud = cyoda.getIceTypePlus(
    yCloud=resCloud,
    nList=nList,
    path="runOne/",
    firstFrame=1,
    isSlice=False,
    outputFileName="chillPlus.txt",
)

# Dump the rescloud which currently has CHILL Plus classifications
cyoda.writeDump(
    yCloud=resCloud,
    path="runOne/",
    outFile="waterChillP.lammpstrj",
)

# Average Q6 (cloud,slice)
avgQ6 = cyoda.getq6(
    yCloud=resCloud,
    nList=nList,
    isSlice=False,
)

# Modification (cloud,q6)
resCloud = cyoda.reclassifyWater(
    yCloud=resCloud,
    q6=avgQ6,
)

# Post reclassification writeOut
cyoda.printIceType(
    yCloud=resCloud,
    path="runOne/",
    firstFrame=1,
    isSlice=False,
    outputFileName="chillPlus.txt",
)

# Dump the rescloud which now has the supaa CHILL Plus Trajectory
cyoda.writeDump(
    yCloud=resCloud,
    path="runOne/",
    outFile="waterSupaaP.lammpstrj",
)

# Get the largest ice cluster. Here, iceNeighbourList is the neighbour list by index.
cyoda.clusterAnalysis(
    path="runOne/",
    iceCloud=solCloud,
    yCloud=resCloud,
    nList=nList,
    iceNeighbourList=iceList,
    cutoff=3.5,
    firstFrame=1,
    bopAnalysis="q6",
)

# Recenter the cluster such that the centroid is at the center of the simulation box
cyoda.recenterClusterCloud(
    iceCloud=solCloud,
    nList=iceList,
)

# Dump the recentered largest ice cluster
cyoda.writeDump(
    yCloud=resCloud,
    path="runOne/",
    outFile="largestIce.lammpstrj",
)
