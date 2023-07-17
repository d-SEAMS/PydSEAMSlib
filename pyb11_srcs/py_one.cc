#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../subprojects/seams-core/src/include/internal/mol_sys.hpp"
#include "../subprojects/seams-core/src/include/internal/seams_input.hpp"

#include <fmt/core.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(cyoda, m) {
    py::class_<molSys::Point<double>>(m, "PointDouble")
        .def(py::init<>())
        .def_readwrite("c_type", &molSys::Point<double>::type)
        .def_readwrite("molID", &molSys::Point<double>::molID)
        .def_readwrite("atomID", &molSys::Point<double>::atomID)
        .def_readwrite("iceType", &molSys::Point<double>::iceType)
        .def_readwrite("x", &molSys::Point<double>::x)
        .def_readwrite("y", &molSys::Point<double>::y)
        .def_readwrite("z", &molSys::Point<double>::z)
        .def_readwrite("c_ij", &molSys::Point<double>::c_ij)
        .def_readwrite("inSlice", &molSys::Point<double>::inSlice)
        .def("__repr__",
             [](const molSys::Point<double> &a) {
                 std::uintptr_t ptr_val = std::uintptr_t(&a);
                 return fmt::format("<PointDouble mem_loc:{:x}>", static_cast<uint>(ptr_val));
             })
        .def("__str__", [](const molSys::Point<double> &a) {
            return fmt::format(
                "x: {} y: {} z: {} type: {} molID: {} atomID: {} inSlice: {}",
                a.x,
                a.y,
                a.z,
                a.type,
                a.molID,
                a.atomID,
                a.inSlice);
        });

    py::enum_<molSys::bond_type>(m, "BondType")
        .value("staggered", molSys::bond_type::staggered)
        .value("eclipsed", molSys::bond_type::eclipsed)
        .value("out_of_range", molSys::bond_type::out_of_range);

    py::enum_<molSys::atom_state_type>(m, "AtomStateType")
        .value("cubic", molSys::atom_state_type ::cubic)
        .value("hexagonal", molSys::atom_state_type ::hexagonal)
        .value("water", molSys::atom_state_type ::water)
        .value("interfacial", molSys::atom_state_type ::interfacial)
        .value("clathrate", molSys::atom_state_type ::clathrate)
        .value("interClathrate", molSys::atom_state_type ::interClathrate)
        .value("unclassified", molSys::atom_state_type ::unclassified)
        .value("reCubic", molSys::atom_state_type ::reCubic)
        .value("reHex", molSys::atom_state_type ::reHex);

    py::class_<molSys::Result>(m, "Result")
        .def(py::init<>())
        .def_readwrite("classifier", &molSys::Result::classifier)
        .def_readwrite("c_value", &molSys::Result::c_value);

    py::class_<molSys::PointCloud<molSys::Point<double>, double>>(m, "PointCloudDouble")
        .def(py::init<>())
        .def_readwrite("pts", &molSys::PointCloud<molSys::Point<double>, double>::pts)
        .def_readwrite("currentFrame",
                       &molSys::PointCloud<molSys::Point<double>, double>::currentFrame)
        .def_readwrite("nop", &molSys::PointCloud<molSys::Point<double>, double>::nop)
        .def_readwrite("box", &molSys::PointCloud<molSys::Point<double>, double>::box)
        .def_readwrite("boxLow", &molSys::PointCloud<molSys::Point<double>, double>::boxLow)
        .def_readwrite("idIndexMap",
                       &molSys::PointCloud<molSys::Point<double>, double>::idIndexMap);

    m.def("readXYZ",
          &sinp::readXYZ,
          "A function to populate a PointCloudDouble with data from a file",
          py::arg("filename")); /* std::string */
    // "yCloud"_a /*  molSys::PointCloud<molSys::Point<double>, double>* */
    // m.def("clearPointCloud", &molSys::clearPointCloud,
    //"A function to populate a PointCloudDouble with data from a file",
    // py::arg("yCloud") /*  molSys::PointCloud<molSys::Point<double>, double>* */
    //);

    m.def("readLammpsTrjreduced",
          &sinp::readLammpsTrjreduced,
          py::arg("filename"),
          "targetFrame"_a,
          "typeI"_a,
          "isSlice"_a,
          "coordLow"_a,
          "coordHigh"_a);

    m.def("readLammpsTrjO",
          &sinp::readLammpsTrjO,
          py::arg("filename"),
          "targetFrame"_a,
          "typeO"_a,
          "isSlice"_a,
          "coordLow"_a,
          "coordHigh"_a);
}
