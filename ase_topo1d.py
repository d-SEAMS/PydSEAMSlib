from pyseams import cyoda
import pprint
from pyseams.adapters import _ase

from ase.io import read as aseread
trajectory ="subprojects/seams-core/input/traj/exampleTraj.lammpstrj"

#Get the frame
resCloud = cyoda.readLammpsTrjreduced(
          filename = trajectory,
          targetFrame = 1,
          typeI = 2, # oxygenAtomType
          isSlice = False,
          coordLow = [0,0,0],
          coordHigh = [0,0,0],
)
# pprint.pprint(dir(resCloud))

atms = aseread(trajectory)
# In ASE, we want to work with atomic numbers instead of LAMMPS types
atms = _ase.map_2({1: 'H', 2: 'O'}, atms)
only_O_mask = [x.symbol == 'O' for x in atms]
pcd = _ase.to_pointcloud(atms[only_O_mask])

# #Calculate the neighborlist by ID
# nList = cyoda.neighListO(
#     rcutoff = 3.5,
#     yCloud = resCloud,
#     typeI = 2, #oxygenAtomType
# )


# print(pprint.pformat(nList))
