#include <pybind11/pybind11.h>
#include "../subprojects/seams-core/src/include/internal/mol_sys.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cyoda, m) {
  py::class_<molSys::Point<double>>(m, "PointDouble")
	      .def(py::init<>())
	      .def_readwrite("c_type", &molSys::Point<double>::type);
}

