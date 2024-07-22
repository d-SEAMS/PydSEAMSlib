from pyseams.cyoda import PointDouble, PointCloudDouble
import numpy as np
import ase


def map_LAMMPS_ID_to_atomic_numbers(dict_map: dict(), atms_a: ase.Atoms):
    pass


def map_1(dict_map: dict(), atms_a: ase.Atoms):
    for key, value in dict_map.items():
        for pt in atms_a:
            if pt.number == key:
                pt.symbol = value
    return atms_a


def map_2(dict_map: dict(), atms_a: ase.Atoms):
    for pt in atms_a:
        if pt.number in dict_map.keys():
            pt.symbol = dict_map[pt.number]
    return atms_a


def to_pointcloud(atms_a: ase.Atoms):
    _pcd = PointCloudDouble()
    # TODO(ruhila): Assumes a rectangular box
    _pcd.box = np.diag(np.asarray(atms_a.get_cell()))
    _pcd.nop = len(atms_a)
    return _pcd
