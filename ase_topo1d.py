from pyseams import cyoda
from pyseams.adapters import _ase
import pprint
import numpy as np

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
# In ASE, we want to work with atomic symbols instead of LAMMPS types
lammps_to_ase = {1: 'H', 2: 'O'}
atms = _ase.map_2(lammps_to_ase, atms)
only_O_mask = [x.symbol == 'O' for x in atms]
# user has to provide proper molID for inslice, each molecule must have one molID. for eg: In H2O, both H atoms and O atom must have same molID. 
# if one wants molHID, they can just change the last ,1) as ,2) in  the following expression:
molOID = np.repeat(np.arange(1,sum(only_O_mask)+1),1) 
pcd = _ase.to_pointcloud(atms, lammps_to_ase, only_O_mask, molOID)

# #Calculate the neighborlist by ID
# nList = cyoda.neighListO(
#     rcutoff = 3.5,
#     yCloud = resCloud,
#     typeI = 2, #oxygenAtomType
# )


# print(pprint.pformat(nList))
