/**
 *  @file   py_one.cc
 *  @author Ruhila
 *  @date   2024-05-25
 *  @brief  Documentation for python bindings (Pyseams)
 *
 *  Created on: 2024-05-25
 *
 *  This is part of the Pyseams project in gsoc24. This documentation helps understanding the python bindings.
 *  The bindings are created using the pybind11 library. 
 *  @ref https://github.com/d-SEAMS/pyseams
 *  @ref https://pyseams.surge.sh
 *  @ref https://pybind11.readthedocs.io 
 *  @ref https://docs.dseams.info
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fmt/core.h>

#include "../subprojects/seams-core/src/include/internal/mol_sys.hpp"
#include "../subprojects/seams-core/src/include/internal/seams_input.hpp"
#include "../subprojects/seams-core/src/include/internal/neighbours.hpp"
#include "../subprojects/seams-core/src/include/internal/bond.hpp"
#include "../subprojects/seams-core/src/include/internal/franzblau.hpp"
#include "../subprojects/seams-core/src/include/internal/ring.hpp"
#include "../subprojects/seams-core/src/include/internal/topo_one_dim.hpp"
#include "../subprojects/seams-core/src/include/internal/bulkTUM.hpp"
#include "../subprojects/seams-core/src/include/internal/selection.hpp"
#include "../subprojects/seams-core/src/include/internal/topo_two_dim.hpp"




namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(cyoda, m) {
 /**
 * @brief This function defines all the point double class.
 *
 */
    py::class_<molSys::Point<double>>(m, "PointDouble") 
        .def(py::init<>())
        .def_readwrite("c_type", &molSys::Point<double>::type) /**< PointDouble value 1 */
        .def_readwrite("molID", &molSys::Point<double>::molID) /**< PointDouble value 2 */
        .def_readwrite("atomID", &molSys::Point<double>::atomID) /**< PointDouble value 3 */
        .def_readwrite("iceType", &molSys::Point<double>::iceType) /**< PointDouble value 4 */
        .def_readwrite("x", &molSys::Point<double>::x) /**< PointDouble value 5 */
        .def_readwrite("y", &molSys::Point<double>::y) /**< PointDouble value 6 */
        .def_readwrite("z", &molSys::Point<double>::z) /**< PointDouble value 7 */
        .def_readwrite("c_ij", &molSys::Point<double>::c_ij) /**< PointDouble value 8 */
        .def_readwrite("inSlice", &molSys::Point<double>::inSlice) /**< PointDouble value 9 */
        .def("__repr__",
             [](const molSys::Point<double> &self_C) {
                 std::uintptr_t ptr_val = std::uintptr_t(&self_C);
                 return fmt::format("<PointDouble mem_loc:{:x}>", static_cast<uint>(ptr_val));
             })
        .def("__str__", [](const molSys::Point<double> &self_C) {
            return fmt::format("x: {} y: {} z: {} type: {} molID: {} atomID: {} inSlice: {}",
                               self_C.x,
                               self_C.y,
                               self_C.z,
                               self_C.type,
                               self_C.molID,
                               self_C.atomID,
                               self_C.inSlice);
        });

/**
 * @brief This function defines all BondType class.
 *
 */
    py::enum_<molSys::bond_type>(m, "BondType")
        .value("staggered", molSys::bond_type::staggered)
        .value("eclipsed", molSys::bond_type::eclipsed)
        .value("out_of_range", molSys::bond_type::out_of_range);

/**
 * @brief This function defines all AtomStateType class.
 *
 */
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

/**
 * @brief This function defines all Result class.
 *
 */

    py::class_<molSys::Result>(m, "Result")
        .def(py::init<>())
        .def_readwrite("classifier", &molSys::Result::classifier)
        .def_readwrite("c_value", &molSys::Result::c_value)
        .def("__repr__",
             [](const molSys::Result &self_C) {
                 std::uintptr_t ptr_val = std::uintptr_t(&self_C);
                 return fmt::format("<Result mem_loc:{:x}>", static_cast<uint>(ptr_val));
             })
        .def("__str__", [](const molSys::Result &self_C) {
            return fmt::format("classifier: {} c_value: {}", self_C.classifier, self_C.c_value);
        });

/**
 * @brief This function defines all PointCloudDouble class.
 *
 */

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

/**
 * @details This is a python binding for the function readXYZ which reads in atom coordinates from an XYZ file.
 *  or reads in an XYZ file
 *  @param[in] filename Filename of the trajectory.
 */

    m.def("readXYZ",
          &sinp::readXYZ,
          "A function to populate a PointCloudDouble with data from a file",
          py::arg("filename"));

/**
 * @details This is a python binding for the function readLammpsTrjreduced which reads in only atoms of the desired type and ignores all atoms which are not in the slice as well.
 * @param[in] filename	The name of the lammps trajectory file to be read in
 * @param[in] targetFrame	The frame number whose information will be read in
 * @param[in] typeI	The type ID of the desired type of atoms
 * @param[in] isSlice	This decides whether a slice will be created or not
 * @param[in] coordLow	Contains the lower limits of the slice, if a slice is to be created
 * @param[in] coordHigh	Contains the upper limits of the slice, if a slice is to be created
 */
    m.def("readLammpsTrjreduced",
          &sinp::readLammpsTrjreduced,
          "A Function that reads in only atoms of the desired type and ignores all atoms which "
          "are not in the slice as well",
          py::arg("filename"),
          "targetFrame"_a,
          "typeI"_a,
          "isSlice"_a,
          "coordLow"_a,
          "coordHigh"_a);

/**
 * @details This is a python binding for the function readLammpsTrjO which reads in atom in a specified frame.
 * @param[in] filename	The name of the lammps trajectory file to be read in
 * @param[in] targetFrame	The frame number whose information will be read in
 * @param[in] typeO	The type ID of the Oxygen atoms
 * @param[in] isSlice	This decides whether a slice will be created or not
 * @param[in] coordLow	Contains the lower limits of the slice, if a slice is to be created
 * @param[in] coordHigh	Contains the upper limits of the slice, if a slice is to be created
 */
    m.def("readLammpsTrjO",
          &sinp::readLammpsTrjO,
          "A Function for reading oxygen atom in a specified frame",
          py::arg("filename"),
          "targetFrame"_a,
          "typeO"_a,
          "isSlice"_a,
          "coordLow"_a,
          "coordHigh"_a);

    m.def("readLammpsTrj",
          &sinp::readLammpsTrj,
          "A Function for reading in a specified frame",
          py::arg("filename"),
          "targetFrame"_a,
          "isSlice"_a,
          "coordLow"_a,
          "coordHigh"_a);

    m.def("readBonds",
        &sinp::readBonds,
        "A function that Reads bonds into a vector of vectors from a file with a specific format",
        py::arg("filename"));

    m.def("atomInSlice",
          &sinp::atomInSlice,
          "If this is 3 then the particle is inside the volume slice",
          py::arg("x"),
          "y"_a,
          "z"_a,
          "coordLow"_a,
          "coordHigh"_a);

    m.def("clearNeighbourList",
        &nneigh::clearNeighbourList,
        "Erases memory for a vector of vectors for the neighbour list",
        py::arg("nList"));

    m.def("getNewNeighbourListByIndex",
        &nneigh::getNewNeighbourListByIndex,
        "Gets a neighbour list by index, according to a pointCloud given as the input",
        py::arg("yCloud"),
        "cutoff"_a);

    m.def("halfNeighList",
        &nneigh::halfNeighList,
        "Inefficient O(n^2) implementation of neighbour lists ",
        py::arg("yCloud"),
        "rcutoff"_a,
        "typeI"_a);

    m.def("neighbourListByIndex",
        &nneigh::neighbourListByIndex,
        "Converts the neighbour list build with atom IDs into a neighbour list of atom indices, according to the pointCloud",
        py::arg("yCloud"),
        "nList"_a);

    m.def("neighList",
          &nneigh::neighList,
          "All these functions use atom IDs and not indices",
          py::arg("yCloud"),
          "rcutoff"_a,
          "typeI"_a,
          "typeJ"_a);

    m.def("neighListO",
        &nneigh::neighListO,
        "Inefficient O(n^2) implementation of neighbour lists ",
        py::arg("yCloud"),
        "rcutoff"_a,
        "typeI"_a);

    m.def("createBondsFromCages",
          &bond::createBondsFromCages,
          "Creates a vector of vectors containing bond connectivity information from the rings vector of vectors and cage information",
          py::arg("rings"),
          "cageList"_a,
          "type"_a,
          "nRings"_a);

    m.def("getHbondDistanceOH",
          &bond::getHbondDistanceOH,
          "Calculates the distance of the hydrogen bond between O and H (of different atoms)",
          py::arg("oCloud"),
          "hCloud"_a,
          "oAtomIndex"_a,
          "hAtomIndex"_a);

/**
 * @details This is the python binding for the function populateHbonds which creates a vector of vectors containing the hydrogen bond connectivity information. 
 * It also decides the existence of the hydrogen bond depending on the O–O and O–H vectors from the neighbour list already constructed.
 *  @param[in] filename Filename of the trajectory, with the hydrogen and oxygen coordinates
 *  @param[in] yCloud The input molSys::PointCloud for the oxygen atoms only
 *  @param[in] nList Row-ordered neighbour list by atom ID
 *  @param[in] targetFrame The target or current frame number (starts from 1) and is not the timestep value
 */
    m.def("populateHbonds",
          &bond::populateHbonds,
          "Create a vector of vectors containing the hydrogen bond connectivity information. Decides the existence of the hydrogen bond depending on the O–O and O–H vectors from the neighbour list already constructed",
          py::arg("yCloud"),
          "filename"_a,
          "targetFrame"_a,
          "Htype"_a,
          "nList"_a);

    m.def("populateHbondsWithInputClouds",
        &bond::populateHbondsWithInputClouds,
        "Create a vector of vectors (similar to the neighbour list conventions) containing the hydrogen bond connectivity information. Decides the existence of the hydrogen bond depending on the O–O and O–H vectors from the neighbour list already constructed, taking a PointCloud for the H atoms as input",
        py::arg("yCloud"),
        "hCloud"_a,
        "nList"_a);

    m.def("trimBonds",
        &bond::trimBonds,
        "Remove duplicate bonds",
        py::arg("bonds"));

    m.def("clearGraph",
        &primitive::clearGraph,
        "Function for clearing vectors in Graph after multiple usage",
        py::arg("currentGraph"));

    m.def("countAllRingsFromIndex",
        &primitive::countAllRingsFromIndex,
        "Creates a vector of vectors of all possible rings",
        py::arg("neighHbondList"),
        "maxDepth"_a); 

    m.def("findRings",
          &primitive::findRings,
          "Main function that searches for all rings",
          py::arg("fullGraph"),
          "v"_a,
          "visited"_a,
          "maxDepth"_a,
          "depth"_a,
          "root"_a);

    m.def("populateGraphFromIndices",
        &primitive::populateGraphFromIndices,
        "Creates a graph object and fills it with the information from a neighbour list of INDICES NOT ATOM IDs created before",
        py::arg("nList"));

    m.def("populateGraphFromNListID",
        &primitive::populateGraphFromNListID,
        "Creates a graph object and fills it with the information from a neighbour list and pointCloud created before",
        py::arg("yCloud"),
        "neighHbondList"_a);

    m.def("removeNonSPrings",
        &primitive::removeNonSPrings,
        "Removes the non-SP rings, using the Franzblau shortest path criterion",
        py::arg("fullGraph"));

    m.def("restoreEdgesFromIndices",
        &primitive::restoreEdgesFromIndices,
        "Re-fills the neighbour lists of a graph object from a neighbour list of INDICES NOT ATOM IDs created before",
        py::arg("fullGraph"),
        "nList"_a);

    m.def("ringNetwork",
        &primitive::ringNetwork,
        "Returns a vector of vectors containing the rings",
        py::arg("maxDepth"),
        "nList"_a);

    m.def("shortestPath",
          &primitive::shortestPath,
          "Calculates the shortest path",
          py::arg("fullGraph"),
          "v"_a,
          "visited"_a,
          "maxDepth"_a,
          "depth"_a,
          "path"_a,
          "goal"_a);

    m.def("assignPolygonType",
        &ring::assignPolygonType,
        "Assign an atomType (equal to the number of nodes in the ring) given n-membered rings",
        py::arg("rings"),
        "atomTypes"_a,
        "nRings"_a);

    m.def("assignPrismType",
          &ring::assignPrismType,
          "Assign an atomType (equal to the number of nodes in the ring) given a vector with a list of indices of rings comprising the prisms",
          py::arg("rings"),
          "listPrism"_a,
          "ringSize"_a,
          "ringType"_a,
          "atomTypes"_a,
          "atomState"_a);

    m.def("atomsFromCages",
        &tum3::atomsFromCages,
        "Gets the atoms in the cages of a given cluster",
        py::arg("rings"),
        "cageList"_a,
        "clusterCages"_a);

    m.def("atomsInSingleSlice",
          &gen::atomsInSingleSlice,
          "Given a pointCloud set the inSlice bool for every atom, if the atoms are inside the specified (single) region",
          py::arg("yCloud"),
          "clearPreviousSliceSelection"_a,
          "coordLow"_a,
          "coordHigh"_a);

    m.def("averageRMSDatom",
        &tum3::averageRMSDatom,
        "Average the RMSD per atom",
        py::arg("rmsdPerAtom"),
        "noOfCommonAtoms"_a);
 
    m.def("basalPrismConditions",
        &ring::basalPrismConditions,
        "Tests whether two rings are basal rings (true) or not (false) for a prism (strict criterion",
        py::arg("nList"),
        "basal1"_a,
        "basal2"_a);

    m.def("buildRefDDC",
        &tum3::buildRefDDC,
        "Build a reference Double-Diamond cage, reading in from a template XYZ file",
        py::arg("fileName"));

    m.def("buildRefHC",
        &tum3::buildRefHC,
        "Build a reference Hexagonal cage, reading in from a template XYZ file",
        py::arg("fileName"));

    m.def("clearRingList",
        &ring::clearRingList,
        "Erases memory for a vector of vectors for a list of rings",
        py::arg("rings"));

    m.def("clusterCages",
          &tum3::clusterCages,
          "Clustering Clusters cages using the Stillinger algorithm and prints out individual XYZ files of clusters.",
          py::arg("yCloud"),
          "path"_a,
          "rings"_a,
          "cageList"_a,
          "numHC"_a,
          "numDDC"_a);
 
    m.def("commonElementsInThreeRings",
        &ring::commonElementsInThreeRings,
        "Common elements in 3 rings",
        py::arg("ring1"),
        "ring2"_a,
        "ring3"_a);

    m.def("compareRings",
        &ring::compareRings,
        "Compares two disordered vectors and checks to see if they contain the same elements",
        py::arg("ring1"),
        "ring2"_a);

    m.def("deformedPrismTypes",
        &ring::deformedPrismTypes,
        "Get the atom type values for deformed prisms",
        py::arg("atomState"),
        "atomTypes"_a,
        "maxDepth"_a);

    m.def("discardExtraTetragonBlocks",
        &ring::discardExtraTetragonBlocks,
        "Checks whether two 4-membered rings are parallel in one dimension or not to prevent overcounting",
        py::arg("basal1"),
        "basal2"_a,
        "yCloud"_a);

    m.def("findPrisms",
        &ring::findPrisms,
        "Find out which rings are prisms. Returns a vector containing all the ring IDs which are prisms",
        py::arg("rings"),
        "ringType"_a,
        "nPerfectPrisms"_a,
        "nImperfectPrisms"_a,
        "nList"_a,
        "rmsdPerAtom"_a,
        "doShapeMatching"_a,
        "yCloud"_a);

    m.def("findsCommonElements",
        &ring::findsCommonElements,
        "Returns the common elements of two rings",
        py::arg("ring1"),
        "ring2"_a);

    m.def("findTripletInRing",
        &ring::findTripletInRing,
        "Searches a particular ring for a triplet",
        py::arg("ring"),
        "triplet"_a);

    m.def("getEdgeMoleculesInRings",
        &ring::getEdgeMoleculesInRings,
        py::arg("rings"),
        "oCloud"_a,
        "yCloud"_a,
        "identicalCloud"_a,
        "coordLow"_a,
        "coordHigh"_a);

    m.def("getPointCloudOneAtomType",
        &gen::getPointCloudOneAtomType,
        "Given a pointCloud containing certain atom types, this returns a pointCloud containing atoms of only the desired type",
        py::arg("yCloud"),
        "outCloud"_a,
        "atomTypeI"_a,
        "isSlice"_a,
        "coordLow"_a,
        "coordHigh"_a);

    m.def("getSingleRingSize",
        &ring::getSingleRingSize,
        "Returns a vector of vectors of rings of a single size.",
        py::arg("rings"),
        "ringSize"_a);

    m.def("hasCommonElements",
        &ring::hasCommonElements,
        "Check to see if two vectors have common elements or not True, if common elements are present and false if there are no common element",
        py::arg("ring1"),
        "ring2"_a);

//   m.def("keepAxialRingsOnly",
//        &ring::keepAxialRingsOnly,
//        "Saves only axial rings out of all possible rings",
//        py::arg("rings"),
//        "yCloud"_a);

    m.def("moleculesInSingleSlice",
        &gen::moleculesInSingleSlice,
        "Given a pointCloud set the inSlice bool for every atom, if the molecules are inside the specified (single) region. If even one atom of a molecule is inside the region, then all atoms of that molecule will be inside the region (irrespective of type)",
        py::arg("yCloud"),
        "clearPreviousSliceSelection"_a,
        "coordLow"_a,
        "coordHigh"_a);
  
    m.def("polygonRingAnalysis",
        &ring::polygonRingAnalysis,
        "Find out which rings are prisms, looping through all ring sizes upto the maxDepth The input ringsAllSizes array has rings of every size",
        py::arg("path"),
        "rings"_a,
        "nList"_a,
        "yCloud"_a,
        "maxDepth"_a,
        "sheetArea"_a,
        "firstFrame"_a);

    m.def("printSliceGetEdgeMoleculesInRings",
        &ring::printSliceGetEdgeMoleculesInRings,
        py::arg("path"),
        "rings"_a,
        "oCloud"_a,
        "yCloud"_a,
        "coordLow"_a,
        "coordHigh"_a,
        "identicalCloud"_a);

    m.def("prismAnalysis",
        &ring::prismAnalysis,
        "Find out which rings are prisms, looping through all ring sizes upto the maxDepth The input ringsAllSizes array has rings of every size",
        py::arg("path"),
        "rings"_a,
        "nList"_a,
        "yCloud"_a,
        "maxDepth"_a,
        "atomID"_a,
        "firstFrame"_a,
        "currentFrame"_a,
        "doShapeMatching"_a);

    m.def("relaxedPrismConditions",
        &ring::relaxedPrismConditions,
        "Two candidate basal rings of a prism block should have at least one bond between them",
        py::arg("nList"),
        "basal1"_a,
        "basal2"_a);

    m.def("rmAxialTranslations",
        &ring::rmAxialTranslations,
        "Shift the entire ice nanotube and remove axial translations",
        py::arg("yCloud"),
        "atomID"_a,
        "firstFrame"_a,
        "currentFrame"_a);

    m.def("setAtomsWithSameMolID",
        &gen::setAtomsWithSameMolID,
        "Given a particular molecule ID and a pointCloud set the inSlice bool for all atoms, with that molecule ID",
        py::arg("yCloud"),
        "molIDAtomIDmap"_a,
        "molID"_a,
        "inSliceValue"_a);
 
    m.def("shapeMatchDDC",
        &tum3::shapeMatchDDC,
        "Shape-matching for a target DDC",
        py::arg("yCloud"),
        "refPoints"_a,
        "cageList"_a,
        "cageIndex"_a,
        "rings"_a,
        "quat"_a,
        "rmsd"_a);

    m.def("shapeMatchHC",
        &tum3::shapeMatchHC,
        "Shape-matching for a target HC",
        py::arg("yCloud"),
        "refPoints"_a,
        "cageUnit"_a,
        "rings"_a,
        "nList"_a,
        "quat"_a,
        "rmsd"_a);

    m.def("topoBulkCriteria",
        &tum3::topoBulkCriteria,
        "Topological network methods Finds the HCs and DDCs for the system",
        py::arg("path"),
        "rings"_a,
        "nList"_a,
        "yCloud"_a,
        "firstFrame"_a,
        "numHC"_a,
        "numDDC"_a,
        "ringType"_a);

    m.def("topoUnitMatchingBulk",
        &tum3::topoUnitMatchingBulk,
        "Topological unit matching for bulk water. If printClusters is true, individual clusters of connected cages are printed",
        py::arg("path"),
        "rings"_a,
        "nList"_a,
        "yCloud"_a,
        "firstFrame"_a,
        "printClusters"_a,
        "onlyTetrahedral"_a);

    m.def("updateRMSDatom",
        &tum3::updateRMSDatom,
        "Calulate the RMSD for each ring, using RMSD values (rmsd) obtained from the shape-matching of each cage",
        py::arg("rings"),
        "cageUnit"_a,
        "rmsd"_a,
        "rmsdPerAtom"_a,
        "noOfCommonAtoms"_a,
        "atomTypes"_a);
}
