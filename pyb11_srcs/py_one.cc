#include <pybind11/pybind11.h>
#include "../subprojects/seams-core/src/include/internal/mol_sys.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cyoda, m) {
  py::class_<molSys::Point<double>>(m, "PointDouble")
	      .def(py::init<>())
	      .def_readwrite("c_type", &molSys::Point<double>::type)
	      .def_readwrite("molID", &molSys::Point<double>::molID)
	      .def_readwrite("atomID", &molSys::Point<double>::atomID) 
	      .def_readwrite("x", &molSys::Point<double>::x)
	      .def_readwrite("y", &molSys::Point<double>::y)
	      .def_readwrite("z", &molSys::Point<double>::z);

py::enum_<molSys::bond_type>(m, "BondType")
    .value("staggered", molSys::bond_type::staggered)
    .value("eclipsed", molSys::bond_type::eclipsed)
    .value("out_of_range", molSys::bond_type::out_of_range);
}
