/*
** This file is part of d-SEAMS (PydSEAMSlib).
**
** SPDX-License-Identifier: MIT
**
** Copyright (c) 2023--present, d-SEAMS core team
** All rights reserved.
**
** Repo:
** https://github.com/d-SEAMS/PydSEAMSlib
*/
#include <pybind11/stl.h>

#include "../subprojects/seams-core/src/include/internal/bond.hpp"
#include "../subprojects/seams-core/src/include/internal/bop.hpp"
#include "../subprojects/seams-core/src/include/internal/bulkTUM.hpp"
#include "../subprojects/seams-core/src/include/internal/cluster.hpp"
#include "../subprojects/seams-core/src/include/internal/franzblau.hpp"
#include "../subprojects/seams-core/src/include/internal/mol_sys.hpp"
#include "../subprojects/seams-core/src/include/internal/neighbours.hpp"
#include "../subprojects/seams-core/src/include/internal/rdf2d.hpp"
#include "../subprojects/seams-core/src/include/internal/ring.hpp"
#include "../subprojects/seams-core/src/include/internal/seams_input.hpp"
#include "../subprojects/seams-core/src/include/internal/seams_output.hpp"
#include "../subprojects/seams-core/src/include/internal/selection.hpp"
#include "../subprojects/seams-core/src/include/internal/topo_one_dim.hpp"
#include "../subprojects/seams-core/src/include/internal/topo_two_dim.hpp"

#include <cstdint>

#include <fmt/core.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(cyoda, m) {
    m.doc() = R"pbdoc(
        PydSEAMSlib bindings
        --------------------

        .. currentmodule:: pydseamslib.cyoda

        .. autosummary::
           :toctree: _generate

           readXYZ
           readLammpsTrjreduced
           readLammpsTrjO
           readLammpsTrj
           readBonds
           atomInSlice
           clearNeighbourList
           getNewNeighbourListByIndex
           halfNeighList
           neighbourListByIndex
           neighList
           neighListO
           createBondsFromCages
           getHbondDistanceOH
           populateHbonds
           populateHbondsWithInputClouds
           trimBonds
           clearGraph
           countAllRingsFromIndex
           findRings
           populateGraphFromIndices
           populateGraphFromNListID
           removeNonSPrings
           restoreEdgesFromIndices
           ringNetwork
           shortestPath
           assignPolygonType
           assignPrismType
           atomsFromCages
           atomsInSingleSlice
           averageRMSDatom
           basalPrismConditions
           buildRefDDC
           buildRefHC
           clearRingList
           clusterCages
           commonElementsInThreeRings
           compareRings
           deformedPrismTypes
           discardExtraTetragonBlocks
           findPrisms
           findsCommonElements
           findTripletInRing
           getEdgeMoleculesInRings
           getPointCloudOneAtomType
           getSingleRingSize
           getSingleRingSize
           hasCommonElements
           moleculesInSingleSlice
           polygonRingAnalysis
           printSliceGetEdgeMoleculesInRings
           prismAnalysis
           relaxedPrismConditions
           rmAxialTranslations
           setAtomsWithSameMolID
           shapeMatchDDC
           shapeMatchHC
           topoBulkCriteria
           topoUnitMatchingBulk
           updateRMSDatom
           rdf2Danalysis_AA
           bulkPolygonRingAnalysis
           clusterAnalysis
           getCorrelPlus
           getIceTypePlus
           writeDump
           getq6
           reclassifyWater
           printIceType
           recenterClusterCloud

    )pbdoc";

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
             [](const molSys::Point<double> &self_C) {
                 std::uintptr_t ptr_val = std::uintptr_t(&self_C);
                 return fmt::format("<PointDouble mem_loc:{:x}>", static_cast<uint64_t>(ptr_val));
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
        .def_readwrite("c_value", &molSys::Result::c_value)
        .def("__repr__", [](const molSys::Result &self_C) {
            std::uintptr_t ptr_val = std::uintptr_t(&self_C);
            return fmt::format("<Result mem_loc:{:x}>", static_cast<uint64_t>(ptr_val));
        });
    //        .def("__str__", [](const molSys::Result &self_C) {
    //            return fmt::format("classifier: {} c_value: {}", self_C.classifier,
    //            self_C.c_value);
    //        });

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
          py::arg("filename"),
          R"pbdoc(

    A function to populate a PointCloudDouble with data from a file

    Parameters
    ----------

    filename : The name of the file that needs to be read

      )pbdoc");

    m.def("readLammpsTrjreduced",
          &sinp::readLammpsTrjreduced,
          py::arg("filename"),
          "targetFrame"_a,
          "typeI"_a,
          "isSlice"_a,
          "coordLow"_a,
          "coordHigh"_a,
          R"pbdoc(

   A Function that reads in only atoms of the desired type and ignores all atoms which are not in the slice as well.

   Parameters
   ----------

   filename : The name of the file that needs to be read
   targetFrame : The frame number whose information will be read in
   typeI : The type ID of the desired type of atoms
   isSlice : This decides whether a slice will be created or not
   coordLow : Contains the lower limits of the slice, if a slice is to be created
   coordHigh : Contains the upper limits of the slice, if a slice is to be created

      )pbdoc");

    m.def("readLammpsTrjO",
          &sinp::readLammpsTrjO,
          py::arg("filename"),
          "targetFrame"_a,
          "typeO"_a,
          "isSlice"_a,
          "coordLow"_a,
          "coordHigh"_a,
          R"pbdoc(

   A Function for reading oxygen atom in a specified frame

   Parameters
   ----------

   filename : The name of the file that needs to be read
   targetFrame : The frame number whose information will be read in
   typeO : 	The type ID of the Oxygen atoms
   isSlice : This decides whether a slice will be created or not
   coordLow : Contains the lower limits of the slice, if a slice is to be created
   coordHigh : Contains the upper limits of the slice, if a slice is to be created

      )pbdoc");

    m.def("readLammpsTrj",
          &sinp::readLammpsTrj,
          py::arg("filename"),
          "targetFrame"_a,
          "isSlice"_a,
          "coordLow"_a,
          "coordHigh"_a,
          R"pbdoc(

   A Function for reading in a specified frame

   Parameters
   ----------

   filename : The name of the file that needs to be read
   targetFrame : The frame number whose information will be read in
   isSlice : This decides whether a slice will be created or not
   coordLow : Contains the lower limits of the slice, if a slice is to be created
   coordHigh : Contains the upper limits of the slice, if a slice is to be created

      )pbdoc");

    m.def("readBonds",
          &sinp::readBonds,
          py::arg("filename"),
          R"pbdoc(

    Reads bonds into a vector of vectors from a file with a specific format

    Parameters
    ----------

    filename : The name of the file that needs to be read

      )pbdoc");

    m.def("atomInSlice",
          &sinp::atomInSlice,
          py::arg("x"),
          "y"_a,
          "z"_a,
          "coordLow"_a,
          "coordHigh"_a,
          R"pbdoc(

   If this is 3 then the particle is inside the volume slice

   Parameters
   ----------

   x : X co-ordinates of the atom
   y : Y co-ordinates of the atom
   z : Z co-ordinates of the atom
   coordLow : Contains the lower limits of the slice, if a slice is to be created
   coordHigh : Contains the upper limits of the slice, if a slice is to be created

      )pbdoc");

    m.def("clearNeighbourList",
          &nneigh::clearNeighbourList,
          "Erases memory for a vector of vectors for the neighbour list",
          py::arg("nList"),
          R"pbdoc(

    Erases memory for a vector of vectors for the neighbour list, Call this before creating the neighbour list for a new frame

    Parameters
    ----------

    nList : Vector of vectors, of the neighbour list to be erased

      )pbdoc");

    m.def("getNewNeighbourListByIndex",
          &nneigh::getNewNeighbourListByIndex,
          "Gets a neighbour list by index, according to a pointCloud given as the input",
          py::arg("yCloud"),
          "cutoff"_a,
          R"pbdoc(

    Gets a neighbour list by index, according to a pointCloud given as the input

    Parameters
    ----------

    yCloud : The input molSys::PointCloud
    cutoff : Distance cutoff, within which two atoms are neighbours

      )pbdoc");

    m.def("halfNeighList",
          &nneigh::halfNeighList,
          "Inefficient O(n^2) implementation of neighbour lists ",
          py::arg("yCloud"),
          "rcutoff"_a,
          "typeI"_a,
          R"pbdoc(

    Inefficient O(n^2) implementation of neighbour lists

    Parameters
    ----------

    yCloud : The input molSys::PointCloud
    cutoff : Distance cutoff, within which two atoms are neighbours
    typeI : Type ID of the i^th particle type

      )pbdoc");

    m.def("neighbourListByIndex",
          &nneigh::neighbourListByIndex,
          "Converts the neighbour list build with atom IDs into a neighbour list of atom indices, "
          "according to the pointCloud",
          py::arg("yCloud"),
          "nList"_a,
          R"pbdoc(

    Converts the neighbour list build with atom IDs into a neighbour list of atom indices, according to the pointCloud

    Parameters
    ----------

    yCloud : The input molSys::PointCloud
    nList : Full neighbour list, by atom ID

      )pbdoc");

    m.def("neighList",
          &nneigh::neighList,
          "All these functions use atom IDs and not indices",
          py::arg("yCloud"),
          "rcutoff"_a,
          "typeI"_a,
          "typeJ"_a,
          R"pbdoc(

    All these functions use atom IDs and not indices

    Parameters
    ----------

    yCloud : The input molSys::PointCloud
    rcutoff : Distance cutoff, within which two atoms are neighbours
    typeI : Type ID of particles of type I
    typeJ : Type ID of particles of type J

      )pbdoc");

    m.def("neighListO",
          &nneigh::neighListO,
          "Inefficient O(n^2) implementation of neighbour lists ",
          py::arg("rcutoff"),
          "yCloud"_a,
          "typeI"_a,
          R"pbdoc(

    Inefficient O(n^2) implementation of neighbour lists

    Parameters
    ----------

    yCloud : The input molSys::PointCloud
    rcutoff : Distance cutoff, within which two atoms are neighbours
    typeI : Type ID of the i^th particle type

      )pbdoc");

    m.def("createBondsFromCages",
          &bond::createBondsFromCages,
          "Creates a vector of vectors containing bond connectivity information from the rings "
          "vector of vectors and cage information",
          py::arg("rings"),
          "cageList"_a,
          "type"_a,
          "nRings"_a,
          R"pbdoc(

    Creates a vector of vectors containing bond connectivity information from the rings vector of vectors and cage information

    Parameters
    ----------

    rings : Row-ordered vector of vectors atom indices of ring information. Each row is a ring, containing the indices of the particles in that ring
    cageList : A vector of cage::Cage containing a list of HCs or DDCs
    type : The type of cage to get bonds for
    nRings : The total number of rings for all the cages, for the particular cage type

      )pbdoc");

    m.def("getHbondDistanceOH",
          &bond::getHbondDistanceOH,
          "Calculates the distance of the hydrogen bond between O and H (of different atoms)",
          py::arg("oCloud"),
          "hCloud"_a,
          "oAtomIndex"_a,
          "hAtomIndex"_a,
          R"pbdoc(

    Calculates the distance of the hydrogen bond between O and H (of different atoms)

    Parameters
    ----------

    oCloud : The molSys::PointCloud for the oxygen atoms only
    hCloud : The molSys::PointCloud for the hydrogen atoms only
    oAtomIndex : The index (in the oCloud) of the oxygen atom
    hAtomIndex : The index (in the hCloud) of the hydrogen atom

      )pbdoc");

    m.def("populateHbonds",
          &bond::populateHbonds,
          "Create a vector of vectors containing the hydrogen bond connectivity information. "
          "Decides the existence of the hydrogen bond depending on the O–O and O–H vectors from "
          "the neighbour list already constructed",
          py::arg("filename"),
          "yCloud"_a,
          "nList"_a,
          "targetFrame"_a,
          "Htype"_a,
          R"pbdoc(

    Create a vector of vectors containing the hydrogen bond connectivity information. Decides the existence of the hydrogen bond depending on the O–O and O–H vectors from the neighbour list already constructed

    Parameters
    ----------

    filename : Filename of the trajectory, with the hydrogen and oxygen coordinates
    yCloud : The input molSys::PointCloud for the oxygen atoms only
    nList : Row-ordered neighbour list by atom ID
    targetFrame : The target or current frame number (starts from 1) and is not the timestep value
    Htype : The type ID of the hydrogen atoms

      )pbdoc");

    m.def("populateHbondsWithInputClouds",
          &bond::populateHbondsWithInputClouds,
          "Create a vector of vectors (similar to the neighbour list conventions) containing the "
          "hydrogen bond connectivity information. Decides the existence of the hydrogen bond "
          "depending on the O–O and O–H vectors from the neighbour list already constructed, "
          "taking a PointCloud for the H atoms as input",
          py::arg("yCloud"),
          "hCloud"_a,
          "nList"_a,
          R"pbdoc(

    Create a vector of vectors (similar to the neighbour list conventions) containing the hydrogen bond connectivity information. Decides the existence of the hydrogen bond depending on the O–O and O–H vectors from the neighbour list already constructed, taking a PointCloud for the H atoms as input

    Parameters
    ----------

    yCloud : The input **molSys::PointCloud** for the hydrogen atoms only, for the entire system (regardless of whether there is a slice or not)
    hCloud : The input **molSys::PointCloud** for the oxygen atoms only
    nList : Row-ordered neighbour list by atom ID

      )pbdoc");

    m.def("trimBonds",
          &bond::trimBonds,
          "Remove duplicate bonds",
          py::arg("bonds"),
          R"pbdoc(

    Remove duplicate bonds

    Parameters
    ----------

    bonds : Row-ordered vector of vectors of the bond matrix Each row is a ring, containing the indices of the particles in that ring

      )pbdoc");

    m.def("clearGraph",
          &primitive::clearGraph,
          "Function for clearing vectors in Graph after multiple usage",
          py::arg("currentGraph"),
          R"pbdoc(

    Function for clearing vectors in Graph after multiple usage

    Parameters
    ----------

    currentGraph : The cleared Graph

      )pbdoc");

    m.def("countAllRingsFromIndex",
          &primitive::countAllRingsFromIndex,
          "Creates a vector of vectors of all possible rings",
          py::arg("neighHbondList"),
          "maxDepth"_a,
          R"pbdoc(

    Creates a vector of vectors of all possible rings

    Parameters
    ----------

    neighHbondList : Row-ordered neighbour list by atom index [not ID]
    maxDepth : The maximum depth upto which rings will be searched. This means that rings larger than maxDepth will not be generated

      )pbdoc");

    m.def("findRings",
          &primitive::findRings,
          "Main function that searches for all rings",
          py::arg("fullGraph"),
          "v"_a,
          "visited"_a,
          "maxDepth"_a,
          "depth"_a,
          "root"_a,
          R"pbdoc(

    Main function that searches for all rings

    Parameters
    ----------

    fullGraph : 	Graph object containing the vertices [and the neighbour lists]. Vertices may be 'removed' from the Graph
    v : The current vertex being visited or checked. It is added to the list of all vertices visited.
    visited : A vector containing a list of the vertices visited for book-keeping. If the current visited vector fulfills the condition for being a ring, it is added to the rings vector of vector in the Graph
    maxDepth : The maximum depth upto which rings will be searched. This means that rings larger than maxDepth will not be generated
    depth : The current depth. When this function is called for the first time from **primitive::countAllRingsFromIndex**, the depth is initialized to 0. When the depth is greater than or equal to maxDepth, the function exits.
    root : The first vertex, from which the current visited vector [candidate ring] is being grown. This is initialized to a dummy value of -1, when it is called from **primitive::countAllRingsFromIndex**

      )pbdoc");

    m.def("populateGraphFromIndices",
          &primitive::populateGraphFromIndices,
          "Creates a graph object and fills it with the information from a neighbour list of "
          "INDICES NOT ATOM IDs created before",
          py::arg("nList"),
          R"pbdoc(

    Creates a graph object and fills it with the information from a neighbour list of INDICES NOT ATOM IDs created before

    Parameters
    ----------

    nList : The row-ordered neighbour list, containing atom indices (according to the input PointCloud)

      )pbdoc");

    m.def("populateGraphFromNListID",
          &primitive::populateGraphFromNListID,
          "Creates a graph object and fills it with the information from a neighbour list and "
          "pointCloud created before",
          py::arg("yCloud"),
          "neighHbondList"_a,
          R"pbdoc(

    Creates a graph object and fills it with the information from a neighbour list and pointCloud created before

    Parameters
    ----------

    yCloud : The input PointCloud.
    neighHbondList : The row-ordered neighbour list, containing atom IDs, and not the atom indices

      )pbdoc");

    m.def("removeNonSPrings",
          &primitive::removeNonSPrings,
          "Removes the non-SP rings, using the Franzblau shortest path criterion",
          py::arg("fullGraph"),
          R"pbdoc(

    Removes the non-SP rings, using the Franzblau shortest path criterion

    Parameters
    ----------

    fullGraph : The Graph object for the current frame. This also contains the rings vector of vectors, which has all possible rings [possibly including non-SP rings]

      )pbdoc");

    m.def("restoreEdgesFromIndices",
          &primitive::restoreEdgesFromIndices,
          "Re-fills the neighbour lists of a graph object from a neighbour list of INDICES NOT "
          "ATOM IDs created before",
          py::arg("fullGraph"),
          "nList"_a,
          R"pbdoc(

    Re-fills the neighbour lists of a graph object from a neighbour list of INDICES NOT ATOM IDs created before

    Parameters
    ----------

    fullGraph : The Graph object for the current frame. The neighbour lists of component Vertex objects may have been depleted
    nList : The row-ordered neighbour list, containing atom indices [according to the input PointCloud]

      )pbdoc");

    m.def("ringNetwork",
          &primitive::ringNetwork,
          "Returns a vector of vectors containing the rings",
          py::arg("nList"),
          "maxDepth"_a,
          R"pbdoc(

    Returns a vector of vectors containing the rings

    Parameters
    ----------

    nList : Row-ordered neighbour list by index [and NOT the atom ID]
    maxDepth : The maximum depth upto which rings will be searched. This means that rings larger than maxDepth in length will not be generated.

      )pbdoc");

    m.def("shortestPath",
          &primitive::shortestPath,
          "Calculates the shortest path",
          py::arg("fullGraph"),
          "v"_a,
          "visited"_a,
          "maxDepth"_a,
          "depth"_a,
          "path"_a,
          "goal"_a,
          R"pbdoc(

    Calculates the shortest path

    Parameters
    ----------

    fullGraph : The Graph object for the current frame.
    v : The current vertex being checked.
    visited : This vector contains the indices of the vertices visited or checked [for book-keeping].
    maxDepth : The maximum depth or maximum length of the rings.
    depth : The current depth. When this function is called from **primitive::removeNonSPrings**, the depth is initialized as 0
    path : The path or length of the visited points [This basically contains all the indices in the visited vector, excluding the current vertex]
    goal : The first element of the candidate ring being checked [the root node]

      )pbdoc");

    m.def("assignPolygonType",
          &ring::assignPolygonType,
          "Assign an atomType (equal to the number of nodes in the ring) given n-membered rings",
          py::arg("rings"),
          "atomTypes"_a,
          "nRings"_a,
          R"pbdoc(

    Assign an atomType (equal to the number of nodes in the ring) given n-membered rings

    Parameters
    ----------

    rings : The vector of vectors containing the primitive rings, of a particular ring size.
    atomTypes : A vector which contains a type for each atom, depending on it's type as classified by the prism identification scheme.
    nRings : Number of rings

      )pbdoc");

    m.def("assignPrismType",
          &ring::assignPrismType,
          "Assign an atomType (equal to the number of nodes in the ring) given a vector with a "
          "list of indices of rings comprising the prisms",
          py::arg("rings"),
          "listPrism"_a,
          "ringSize"_a,
          "ringType"_a,
          "atomTypes"_a,
          "atomState"_a,
          R"pbdoc(

    Assign an atomType (equal to the number of nodes in the ring) given a vector with a list of indices of rings comprising the prisms

    Parameters
    ----------

    rings : The vector of vectors containing the primitive rings, of a particular ring size
    listPrism : The list of prism blocks found
    ringSize : The current ring size or number of nodes in each ring.
    ringType : A vector which contains a type for each atom, depending on it's type as classified by the prism identification scheme
    atomTypes : A vector which contains a type for each atom, depending on it's type as classified by the prism identification scheme.
    atomState : The state of the atom

      )pbdoc");

    // TODO: read about cagelist and clustercages function
    m.def("atomsFromCages",
          &tum3::atomsFromCages,
          "Gets the atoms in the cages of a given cluster",
          py::arg("rings"),
          "cageList"_a,
          "clusterCages"_a,
          R"pbdoc(

    Gets the atoms in the cages of a given cluster

    Parameters
    ----------

    rings : The vector of vectors containing the primitive rings, of a particular ring size
    cageList : A vector of cage::Cage containing a list of HCs or DDCs
    clusterCages : The cluster of the cages

      )pbdoc");

    m.def("atomsInSingleSlice",
          &gen::atomsInSingleSlice,
          "Given a pointCloud set the inSlice bool for every atom, if the atoms are inside the "
          "specified (single) region",
          py::arg("yCloud"),
          "clearPreviousSliceSelection"_a,
          "coordLow"_a,
          "coordHigh"_a,
          R"pbdoc(

    Given a pointCloud set the inSlice bool for every atom, if the atoms are inside the specified (single) region

    Parameters
    ----------

    yCloud : The given input PointCloud
    clearPreviousSliceSelection : sets all inSlice bool values to false before adding Points to the slice
    coordLow : Contains the lower limits of the slice, if a slice is to be created
    coordHigh : Contains the upper limits of the slice, if a slice is to be created

      )pbdoc");

    // TODO: read about rmsdPerAtom and noOfCommonAtoms function

    m.def("averageRMSDatom",
          &tum3::averageRMSDatom,
          "Average the RMSD per atom",
          py::arg("rmsdPerAtom"),
          "noOfCommonAtoms"_a,
          R"pbdoc(

    Average the RMSD per atom

    Parameters
    ----------

    rmsdPerAtom : The root mean square distance per atom
    noOfCommonAtoms : The number od common atom

      )pbdoc");

    m.def("basalPrismConditions",
          &ring::basalPrismConditions,
          "Tests whether two rings are basal rings (true) or not (false) for a prism (strict "
          "criterion",
          py::arg("nList"),
          "basal1"_a,
          "basal2"_a,
          R"pbdoc(

    Tests whether two rings are basal rings (true) or not (false) for a prism (strict criterion)

    Parameters
    ----------

    nList : Row-ordered neighbour list by atom index
    basal1 : The vector for one of the basal rings
    basal2 : The vector for the other basal ring.

      )pbdoc");

    m.def("buildRefDDC",
          &tum3::buildRefDDC,
          "Build a reference Double-Diamond cage, reading in from a template XYZ file",
          py::arg("fileName"),
          R"pbdoc(

    Build a reference Double-Diamond cage, reading in from a template XYZ file

    Parameters
    ----------

    filename : The name of the file that needs to be read

      )pbdoc");

    m.def("buildRefHC",
          &tum3::buildRefHC,
          "Build a reference Hexagonal cage, reading in from a template XYZ file",
          py::arg("fileName"),
          R"pbdoc(

    Build a reference Hexagonal cage, reading in from a template XYZ file

    Parameters
    ----------

    filename : The name of the file that needs to be read

      )pbdoc");

    m.def("clearRingList",
          &ring::clearRingList,
          "Erases memory for a vector of vectors for a list of rings",
          py::arg("rings"),
          R"pbdoc(

    Erases memory for a vector of vectors for a list of rings

    Parameters
    ----------

    rings : The vector of vectors to be cleared

      )pbdoc");

    // TODO: read about path, numHC, numDDC and cageList function

    m.def("clusterCages",
          &tum3::clusterCages,
          "Clustering Clusters cages using the Stillinger algorithm and prints out individual XYZ "
          "files of clusters.",
          py::arg("yCloud"),
          "path"_a,
          "rings"_a,
          "cageList"_a,
          "numHC"_a,
          "numDDC"_a,
          R"pbdoc(

    Clustering Clusters cages using the Stillinger algorithm and prints out individual XYZ files of clusters.

    Parameters
    ----------

    yCloud : The given input PointCloud
    path : The string to the output directory, in which files will be written out
    rings : The vector of vectors containing the primitive rings, of a particular ring size.
    cageList : A vector of cage::Cage containing a list of HCs or DDCs
    numHC : The number of Hexagonal cage
    numDDC : The number of Double-Diamond cage

      )pbdoc");

    m.def("commonElementsInThreeRings",
          &ring::commonElementsInThreeRings,
          "Common elements in 3 rings",
          py::arg("ring1"),
          "ring2"_a,
          "ring3"_a,
          R"pbdoc(

    Common elements in 3 rings

    Parameters
    ----------

    ring1 : The first ring
    ring2 : The second ring
    ring3 : The third ring

      )pbdoc");

    m.def("compareRings",
          &ring::compareRings,
          "Compares two disordered vectors and checks to see if they contain the same elements",
          py::arg("ring1"),
          "ring2"_a,
          R"pbdoc(

    Compares two disordered vectors and checks to see if they contain the same elements

    Parameters
    ----------

    ring1 : The first ring
    ring2 : The second ring

      )pbdoc");

    m.def("deformedPrismTypes",
          &ring::deformedPrismTypes,
          "Get the atom type values for deformed prisms",
          py::arg("atomState"),
          "atomTypes"_a,
          "maxDepth"_a,
          R"pbdoc(

    Get the atom type values for deformed prisms

    Parameters
    ----------

    atomState : The state of the atom
    atomTypes : A vector which contains a type for each atom, depending on it's type as classified by the prism identification scheme
    maxDepth : The maximum depth or maximum length of the rings.

      )pbdoc");

    m.def("discardExtraTetragonBlocks",
          &ring::discardExtraTetragonBlocks,
          "Checks whether two 4-membered rings are parallel in one dimension or not to prevent "
          "overcounting",
          py::arg("basal1"),
          "basal2"_a,
          "yCloud"_a,
          R"pbdoc(

    Checks whether two 4-membered rings are parallel in one dimension or not to prevent overcounting

    Parameters
    ----------

    basal1 : The first basal ring.
    basal2 : The other candidate basal ring.
    yCloud : The input PointCloud.

      )pbdoc");

    m.def("findPrisms",
          &ring::findPrisms,
          "Find out which rings are prisms. Returns a vector containing all the ring IDs which "
          "are prisms",
          py::arg("rings"),
          "ringType"_a,
          "nPerfectPrisms"_a,
          "nImperfectPrisms"_a,
          "nList"_a,
          "rmsdPerAtom"_a,
          "doShapeMatching"_a,
          "yCloud"_a,
          R"pbdoc(

    Find out which rings are prisms. Returns a vector containing all the ring IDs which are prisms

    Parameters
    ----------

    rings : The input vector of vectors containing the primitive rings of a single ring size [number of nodes]
    ringType : A vector containing a **ring::strucType** value [a classification type] for each ring.
    nPerfectPrisms : The number of perfect prism blocks identified for the number of nodes
    nImperfectPrisms : The number of imperfect prism blocks identified for the number of nodes
    nList : The row-ordered neighbour list [by atom index]
    rmsdPerAtom : The root mean square distance per atom
    doShapeMatching : If one wants to do shape matching
    yCloud : The input PointCloud

      )pbdoc");

    m.def("findsCommonElements",
          &ring::findsCommonElements,
          "Returns the common elements of two rings",
          py::arg("ring1"),
          "ring2"_a,
          R"pbdoc(

    Returns the common elements of two rings

    Parameters
    ----------

    ring1 : The first ring
    ring2 : The second ring

      )pbdoc");

    m.def("findTripletInRing",
          &ring::findTripletInRing,
          "Searches a particular ring for a triplet",
          py::arg("ring"),
          "triplet"_a,
          R"pbdoc(

    Searches a particular ring for a triplet

    Parameters
    ----------

    ring1 : The input ring containing the indices of atoms.
    ring2 : Vector containing the triplet, for whose presence the input ring vector will be checked.


      )pbdoc");

    m.def("getEdgeMoleculesInRings",
          &ring::getEdgeMoleculesInRings,
          py::arg("rings"),
          "oCloud"_a,
          "yCloud"_a,
          "identicalCloud"_a,
          "coordLow"_a,
          "coordHigh"_a,
          R"pbdoc(

    Function that loops through the PointCloud used to construct the neighbour list (used to generate primitive rings) and sets the inSlice bool values of edge atoms which belong to rings that are formed by atoms in the slice.

    Parameters
    ----------

    rings : Vector of vectors of the primitive rings [by index] according to oCloud
    oCloud : PointCloud of O atoms, used to construct the rings vector of vectors
    yCloud : The output PointCloud [may contain more than just the O atoms]
    identicalCloud :  if this is true then oCloud and yCloud are the same
    coordLow : Contains the lower limits of the slice [needed to find all the molecules of oCloud in the slice]
    coordHigh : Contains the upper limits of the slice


      )pbdoc");

    m.def("getPointCloudOneAtomType",
          &gen::getPointCloudOneAtomType,
          "Given a pointCloud containing certain atom types, this returns a pointCloud containing "
          "atoms of only the desired type",
          py::arg("yCloud"),
          "outCloud"_a,
          "atomTypeI"_a,
          "isSlice"_a,
          "coordLow"_a,
          "coordHigh"_a,
          R"pbdoc(

    Given a pointCloud containing certain atom types, this returns a pointCloud containing atoms of only the desired type

    Parameters
    ----------

    atomTypeI : The type ID of the atoms to save in the output PointCloud
    outCloud : The output PointCloud
    yCloud : The given input PointCloud
    isSlice :  This decides whether a slice will be used or not
    coordLow : Contains the lower limits of the slice, if a slice is to be created
    coordHigh : Contains the upper limits of the slice, if a slice is to be created


      )pbdoc");

    m.def("getSingleRingSize",
          &ring::getSingleRingSize,
          "Returns a vector of vectors of rings of a single size.",
          py::arg("rings"),
          "ringSize"_a,
          R"pbdoc(

    Returns a vector of vectors of rings of a single size.

    Parameters
    ----------

    rings : The vector of vectors containing the primitive rings of all sizes
    ringSize : The desired ring size or number of nodes in each ring to be saved.


      )pbdoc");

    m.def("hasCommonElements",
          &ring::hasCommonElements,
          "Check to see if two vectors have common elements or not True, if common elements are "
          "present and false if there are no common element",
          py::arg("ring1"),
          "ring2"_a,
          R"pbdoc(

    Returns a vector of vectors of rings of a single size.

    Parameters
    ----------

    ring1 : The vector of the first ring.
    ring2 : The vector of the second ring.

      )pbdoc");

    //   m.def("keepAxialRingsOnly",
    //        &ring::keepAxialRingsOnly,
    //        "Saves only axial rings out of all possible rings",
    //        py::arg("rings"),
    //        "yCloud"_a);

    m.def("moleculesInSingleSlice",
          &gen::moleculesInSingleSlice,
          "Given a pointCloud set the inSlice bool for every atom, if the molecules are inside "
          "the specified (single) region. If even one atom of a molecule is inside the region, "
          "then all atoms of that molecule will be inside the region (irrespective of type)",
          py::arg("yCloud"),
          "clearPreviousSliceSelection"_a,
          "coordLow"_a,
          "coordHigh"_a,
          R"pbdoc(

    Given a pointCloud set the inSlice bool for every atom, if the atoms are inside the specified (single) region

    Parameters
    ----------

    yCloud : The given input PointCloud
    clearPreviousSliceSelection : sets all inSlice bool values to false before adding Points to the slice
    coordLow : Contains the lower limits of the slice, if a slice is to be created
    coordHigh : Contains the upper limits of the slice, if a slice is to be created

      )pbdoc");

    m.def("polygonRingAnalysis",
          &ring::polygonRingAnalysis,
          "Find out which rings are prisms, looping through all ring sizes upto the maxDepth The "
          "input ringsAllSizes array has rings of every size",
          py::arg("path"),
          "rings"_a,
          "nList"_a,
          "yCloud"_a,
          "maxDepth"_a,
          "sheetArea"_a,
          "firstFrame"_a,
          R"pbdoc(

    Find out which rings are prisms, looping through all ring sizes upto the maxDepth The input ringsAllSizes array has rings of every size

    Parameters
    ----------

    path : The string to the output directory, in which files will be written out.
    rings : Row-ordered vector of vectors for rings of a single type.
    nList : Row-ordered neighbour list by index
    yCloud : The input PointCloud
    maxDepth : The maximum possible size of the primitive rings.
    sheetArea : Area calculated using the two significant dimensions of the quasi-two-dimensional sheet.
    firstFrame : The first frame to be analyzed

      )pbdoc");

    m.def("printSliceGetEdgeMoleculesInRings",
          &ring::printSliceGetEdgeMoleculesInRings,
          py::arg("path"),
          "rings"_a,
          "oCloud"_a,
          "yCloud"_a,
          "coordLow"_a,
          "coordHigh"_a,
          "identicalCloud"_a,
          R"pbdoc(

    Function that loops through the PointCloud used to construct the neighbour list (used to generate primitive rings) and sets the inSlice bool values of edge atoms which belong to rings that are formed by atoms in the slice.

    Parameters
    ----------

    rings : Vector of vectors of the primitive rings [by index] according to oCloud
    path : The string to the output directory, in which files will be written out
    oCloud : PointCloud of O atoms, used to construct the rings vector of vectors
    yCloud : The output PointCloud [may contain more than just the O atoms]
    identicalCloud :  if this is true then oCloud and yCloud are the same
    coordLow : Contains the lower limits of the slice [needed to find all the molecules of oCloud in the slice]
    coordHigh : Contains the upper limits of the slice


      )pbdoc");

    m.def("prismAnalysis",
          &ring::prismAnalysis,
          "Find out which rings are prisms, looping through all ring sizes upto the maxDepth The "
          "input ringsAllSizes array has rings of every size",
          py::arg("path"),
          "rings"_a,
          "nList"_a,
          "yCloud"_a,
          "maxDepth"_a,
          "atomID"_a,
          "firstFrame"_a,
          "currentFrame"_a,
          "doShapeMatching"_a,
          R"pbdoc(

    Find out which rings are prisms, looping through all ring sizes upto the maxDepth The input ringsAllSizes array has rings of every size

    Parameters
    ----------

    rings : Row-ordered vector of vectors for rings of a single type
    path : The string to the output directory, in which files will be written out
    nList : Row-ordered neighbour list by index
    yCloud : The input PointCloud.
    maxDepth : The maximum possible size of the primitive rings
    atomID : The ID of the atom
    firstFrame : The first frame
    currentFrame : The current frame
    doShapeMatching : If one wants to do shape matching


      )pbdoc");

    m.def("relaxedPrismConditions",
          &ring::relaxedPrismConditions,
          "Two candidate basal rings of a prism block should have at least one bond between them",
          py::arg("nList"),
          "basal1"_a,
          "basal2"_a,
          R"pbdoc(

    Two candidate basal rings of a prism block should have at least one bond between them

    Parameters
    ----------

    nList : Row-ordered neighbour list by atom index
    basal1 : The vector for one of the basal rings
    basal2 : The vector for the other basal ring.

      )pbdoc");

    m.def("rmAxialTranslations",
          &ring::rmAxialTranslations,
          "Shift the entire ice nanotube and remove axial translations",
          py::arg("yCloud"),
          "atomID"_a,
          "firstFrame"_a,
          "currentFrame"_a,
          R"pbdoc(

    Shift the entire ice nanotube and remove axial translations

    Parameters
    ----------

    yCloud : The input PointCloud
    atomID : The ID of the atom
    firstFrame : The first frame to be analyzed
    currentFrame : The current frame to be analyzed

      )pbdoc");

    m.def("setAtomsWithSameMolID",
          &gen::setAtomsWithSameMolID,
          "Given a particular molecule ID and a pointCloud set the inSlice bool for all atoms, "
          "with that molecule ID",
          py::arg("yCloud"),
          "molIDAtomIDmap"_a,
          "molID"_a,
          "inSliceValue"_a,
          R"pbdoc(

    Function that loops through a given input PointCloud and sets the inSlice bool for every Point according to whether the molecule is in the specified (single) slice or not. If even one atom of a molecule is inside the region, then all atoms belonging to that molecule should be inside the slice as well (therefore, inSlice would be set to true)

    Parameters
    ----------

    yCloud : The input PointCloud
    molID : The ID of the molecule
    molIDAtomIDmap : Mapping of molecule ID and atom ID
    inSliceValue : It's a bool value

      )pbdoc");

    m.def("shapeMatchDDC",
          &tum3::shapeMatchDDC,
          "Shape-matching for a target DDC",
          py::arg("yCloud"),
          "refPoints"_a,
          "cageList"_a,
          "cageIndex"_a,
          "rings"_a,
          "quat"_a,
          "rmsd"_a,
          R"pbdoc(

    Shape-matching for a target DDC

    Parameters
    ----------

    yCloud : The input PointCloud
    refPoints : The point set of the reference system (or right system). This is a (n*3) Eigen matrix. Here, n is the number of particles.
    cageList : A vector of cage::Cage containing a list of HCs or DDCs
    cageIndex : A vector of cage::Cage containing a Index of HCs or DDCs
    rings : Row-ordered vector of vectors for rings of a single type
    quat : The quaternion for the optimum rotation of the left system into the right system.
    rmsd : Root mean square distance

      )pbdoc");

    m.def("shapeMatchHC",
          &tum3::shapeMatchHC,
          "Shape-matching for a target HC",
          py::arg("yCloud"),
          "refPoints"_a,
          "cageUnit"_a,
          "rings"_a,
          "nList"_a,
          "quat"_a,
          "rmsd"_a,
          R"pbdoc(

    Shape-matching for a target HC

    Parameters
    ----------

    yCloud : The input PointCloud
    refPoints : The point set of the reference system (or right system). This is a (n*3) Eigen matrix. Here, n is the number of particles.
    cageUnit : A vector of cage::Cage containing a units of HCs or DDCs
    nList : Row-ordered neighbour list, by index
    rings : Row-ordered vector of vectors for rings of a single type
    quat : The quaternion for the optimum rotation of the left system into the right system.
    rmsd : Root mean square distance

      )pbdoc");

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
          "ringType"_a,
          R"pbdoc(

    Topological network methods Finds the HCs and DDCs for the system

    Parameters
    ----------

    path : The file path of the output directory to which output files will be written
    rings : Vector of vectors containing the primitive rings. This should contain rings of only size 6
    nList : Row-ordered neighbour list, by index
    yCloud : The input PointCloud, with respect to which the indices in the rings and nList vector of vectors have been saved.
    firstFrame : First frame to be analyzed
    numHC : The number of hexagonal cages
    numDDC : The number of double-diamond cages
    ringType : A vector containing a **ring::strucType** value [a classification type] for each ring.

      )pbdoc");

    m.def("topoUnitMatchingBulk",
          &tum3::topoUnitMatchingBulk,
          "Topological unit matching for bulk water. If printClusters is true, individual "
          "clusters of connected cages are printed",
          py::arg("path"),
          "rings"_a,
          "nList"_a,
          "yCloud"_a,
          "firstFrame"_a,
          "printClusters"_a,
          "onlyTetrahedral"_a,
          R"pbdoc(

    Topological unit matching for bulk water. If printClusters is true, individual clusters of connected cages are printed

    Parameters
    ----------

    path : The file path of the output directory to which output files will be written
    rings : Vector of vectors containing the primitive rings. This should contain rings of only size 6
    nList : Row-ordered neighbour list, by index
    yCloud : The input PointCloud, with respect to which the indices in the rings and nList vector of vectors have been saved.
    firstFrame : First frame to be analyzed
    printClusters : The number of hexagonal cages
    onlyTetrahedral : Flag for only finding DDCs and HCs (true) or also finding PNCs (false)

      )pbdoc");

    m.def("updateRMSDatom",
          &tum3::updateRMSDatom,
          "Calulate the RMSD for each ring, using RMSD values (rmsd) obtained from the "
          "shape-matching of each cage",
          py::arg("rings"),
          "cageUnit"_a,
          "rmsd"_a,
          "rmsdPerAtom"_a,
          "noOfCommonAtoms"_a,
          "atomTypes"_a,
          R"pbdoc(

    Calulate the RMSD for each ring, using RMSD values (rmsd) obtained from the shape-matching of each cage

    Parameters
    ----------

    rings : Vector of vectors containing the primitive rings. This should contain rings of only size 6
    cageUnit : A vector of cage::Cage containing a units of HCs or DDCs
    rmsd : Root mean square distance
    rmsdPerAtom : Root mean square distance per atom
    noOfCommonAtoms : The number of common atoms
    atomTypes : A vector which contains a type for each atom, depending on it's type as classified by the prism identification scheme

      )pbdoc");

    m.def("rdf2Danalysis_AA",
          &rdf2::rdf2Danalysis_AA,
          "calls other functions for initializing, sampling and normalizing the RDF",
          py::arg("path"),
          "rdfValues"_a,
          "yCloud"_a,
          "cutoff"_a,
          "binwidth"_a,
          "firstFrame"_a,
          "finalFrame"_a,
          R"pbdoc(

    Calculates the in-plane RDF for quasi-two-dimensional water, when both the atoms are of the same type

    Parameters
    ----------

    path : The file path of the output directory to which output files will be written
    rdfValues : Vector containing the RDF values
    yCloud : The input PointCloud
    cutoff : Cutoff for the RDF. This should not be greater than half the box length
    binwidth : Width of the bin for histogramming
    firstFrame : The first frame for RDF binning
    finalFrame : The final frame for RDF binning

      )pbdoc");

    m.def("bulkPolygonRingAnalysis",
          &ring::bulkPolygonRingAnalysis,
          "Find out rings in the bulk, looping through all ring sizes upto the maxDepth",
          py::arg("path"),
          "rings"_a,
          "nList"_a,
          "yCloud"_a,
          "maxDepth"_a,
          "firstFrame"_a,
          R"pbdoc(

    Find out rings in the bulk, looping through all ring sizes upto the maxDepth

    Parameters
    ----------

    path : The file path of the output directory to which output files will be written
    rings : Row-ordered vector of vectors for rings of a single type
    nList : Row-ordered neighbour list by index
    yCloud : The input PointCloud
    maxDepth : The maximum possible size of the primitive rings
    firstFrame : The first frame to be analyzed

      )pbdoc");

    m.def("clusterAnalysis",
          &clump::clusterAnalysis,
          "Does the cluster analysis of ice particles in the system.",
          py::arg("path"),
          "iceCloud"_a,
          "yCloud"_a,
          "nList"_a,
          "iceNeighbourList"_a,
          "cutoff"_a,
          "firstFrame"_a,
          "bopAnalysis"_a,
          R"pbdoc(

    Does the cluster analysis of ice particles in the system.

    Parameters
    ----------

    path : The file path of the output directory to which output files will be written
    iceCloud : The **molSys::PointCloud** for the largest ice cluster of the ice-like molecules
    nList : Row-ordered neighbour list by atom ID
    yCloud : The **molSys::PointCloud** for all the particles in the frame, regardless of ice type
    iceNeighbourList : Row-ordered neighbour list by atom index, not ID, according to the iceCloud atoms
    cutoff : Cutoff for the nearest neighbours
    firstFrame : The first frame to be analyzed
    bopAnalysis : This determines which method to use for determining the ice-like nature of the particles. This can be "q6" or "chill", for using the "q6" parameter or CHILL algorithm, respectively

      )pbdoc");

    m.def("getCorrelPlus",
          &chill::getCorrelPlus,
          "Gets c_ij and then classifies bond types according to the CHILL+ algorithm.",
          py::arg("yCloud"),
          "nList"_a,
          "isSlice"_a,
          R"pbdoc(

    Gets c_ij and then classifies bond types according to the CHILL+ algorithm.

    Parameters
    ----------

    yCloud : The output **molSys::PointCloud**
    nList : Row-ordered neighbour list by atom ID
    isSlice : This decides whether there is a slice or not

      )pbdoc");

    m.def("getIceTypePlus",
          &chill::getIceTypePlus,
          "Gets c_ij and then classifies bond types according to the CHILL+ algorithm.",
          py::arg("yCloud"),
          "nList"_a,
          "path"_a,
          "firstFrame"_a,
          "isSlice"_a,
          "outputFileName"_a,
          R"pbdoc(

    Classifies each atom according to the CHILL+ algorithm

    Parameters
    ----------

    yCloud : The output **molSys::PointCloud**
    nList : Row-ordered neighbour list by atom ID
    path : Path to the output directory to which ice types are written out to
    firstFrame : The first frame to be analyzed
    isSlice : This decides whether there is a slice or not
    outputFileName : Name of the output file, to which the ice types will be written out. The default file name is "chillPlus.txt"

      )pbdoc");

    m.def("writeDump",
          &sout::writeDump,
          "Generic function for writing out to a dump file.",
          py::arg("yCloud"),
          "path"_a,
          "outFile"_a,
          R"pbdoc(

    Generic function for writing out to a dump file.

    Parameters
    ----------

    yCloud : The output **molSys::PointCloud**
    path : Path to the output directory to which the info in PairCorrel struct are printed out.
    outFile : The output file to which the info in PairCorrel struct are printed out

      )pbdoc");

    m.def("getq6",
          &chill::getq6,
          "q6 can distinguish between water and ice. Use this for the largest ice cluster.",
          py::arg("yCloud"),
          "nList"_a,
          "isSlice"_a,
          R"pbdoc(

    q6 can distinguish between water and ice. Use this for the largest ice cluster.

    Parameters
    ----------

    yCloud : The output **molSys::PointCloud**
    nList : Row-ordered neighbour list by atom ID
    isSlice : This decides whether there is a slice or not


      )pbdoc");

    m.def("reclassifyWater",
          &chill::reclassifyWater,
          py::arg("yCloud"),
          "q6"_a,
          R"pbdoc(

    Reclassifies atoms which may have been mis-classified as water using the averaged q6 and q3 parameters.

    Parameters
    ----------

    yCloud : The output **molSys::PointCloud**
    q6 : Vector containing the previously calculated averaged q6 values using **chill::getq6**


      )pbdoc");

    m.def("printIceType",
          &chill::printIceType,
          py::arg("yCloud"),
          "path"_a,
          "firstFrame"_a,
          "isSlice"_a,
          "outputFileName"_a,
          R"pbdoc(

    Prints out the iceType for a particular frame onto the terminal. Prints out the **molSys::atom_state_type** per-particle ice type, for a particular frame, to a file.

    Parameters
    ----------

    yCloud : The output **molSys::PointCloud** for the current frame
    path : Path to the output directory to which ice types are written out to
    firstFrame : First frame to be analyzed
    isSlice : Determines whether there is a slice or not
    outputFileName : File name of the output file, to which the per-particle ice types will be written out. The default file name is "superChill.txt"


      )pbdoc");

    m.def("recenterClusterCloud",
          &clump::recenterClusterCloud,
          py::arg("iceCloud"),
          "nList"_a,
          R"pbdoc(

    Recenters the largest ice cluster, by applying a transformation on the largest ice cluster coordinates. Requires the neighbour list BY INDEX

    Parameters
    ----------

    iceCloud : The output **molSys::PointCloud** for the largest ice cluster of the ice-like molecules.
    nList : Row-ordered neighbour list by atom index, for the molSys::PointCloud iceCloud


      )pbdoc");
}
