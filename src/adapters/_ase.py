from pyseams.cyoda import PointDouble, PointCloudDouble
import numpy as np
import ase


def to_pointcloud(atms_a: ase.Atoms):
    _pcd = PointCloudDouble()
    # TODO(ruhila): Assumes a rectangular box
    _pcd.box = np.diag(np.asarray(atms_a.get_cell()))
    return _pcd
