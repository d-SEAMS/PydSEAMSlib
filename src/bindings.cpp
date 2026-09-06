/*
** This file is part of d-SEAMS (PydSEAMSlib).
**
** SPDX-License-Identifier: MIT
**
** Copyright (c) 2023--present, d-SEAMS core team
** All rights reserved.
*/
#include <bond.hpp>
#include <bop.hpp>
#include <bulkTUM.hpp>
#include <cage_affiliation.hpp>
#include <cage_enum.hpp>
#include <cluster.hpp>
#include <density.hpp>
#include <franzblau.hpp>
#include <ira_sofi.hpp>

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/vector.h>
#ifdef SEAMS_HAS_IRA
#    include <Eigen/Core>
#endif
#include <cstdint>
#include <format>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <optional>
#include <rdf.hpp>
#include <rdf2d.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <selection.hpp>
#include <site.hpp>
#include <structure_desc.hpp>
#include <topo_fingerprint.hpp>
#include <topo_one_dim.hpp>
#include <topo_two_dim.hpp>
#include <voronoi_qlm.hpp>

namespace nb = nanobind;

NB_MODULE(yoda, m) {
    m.doc() = "d-SEAMS compiled surface (yoda)";

    nb::class_<molSys::Point<double>>(
        m, "PointDouble", "Per-particle data: coordinates, type, molecule ID, ice classification.")
        .def(nb::init<>())
        .def_rw("c_type", &molSys::Point<double>::type)
        .def_rw("molID", &molSys::Point<double>::molID)
        .def_rw("atomID", &molSys::Point<double>::atomID)
        .def_rw("iceType", &molSys::Point<double>::iceType)
        .def_rw("x", &molSys::Point<double>::x)
        .def_rw("y", &molSys::Point<double>::y)
        .def_rw("z", &molSys::Point<double>::z)
        .def_rw("c_ij", &molSys::Point<double>::c_ij)
        .def_rw("inSlice", &molSys::Point<double>::inSlice)
        .def("__repr__",
             [](const molSys::Point<double> &self_C) {
                 std::uintptr_t ptr_val = std::uintptr_t(&self_C);
                 return std::format("<PointDouble mem_loc:{:x}>", static_cast<uint64_t>(ptr_val));
             })
        .def("__str__", [](const molSys::Point<double> &self_C) {
            return std::format("x: {} y: {} z: {} type: {} molID: {} atomID: {} inSlice: {}",
                               self_C.x,
                               self_C.y,
                               self_C.z,
                               self_C.type,
                               self_C.molID,
                               self_C.atomID,
                               self_C.inSlice);
        });

    nb::enum_<molSys::bond_type>(
        m, "BondType", "Bond classification: staggered, eclipsed, or out_of_range.")
        .value("staggered", molSys::bond_type::staggered)
        .value("eclipsed", molSys::bond_type::eclipsed)
        .value("out_of_range", molSys::bond_type::out_of_range);

    nb::enum_<molSys::atom_state_type>(
        m, "AtomStateType", "Per-atom ice phase classification from CHILL/CHILL+/q6.")
        .value("cubic", molSys::atom_state_type::cubic)
        .value("hexagonal", molSys::atom_state_type::hexagonal)
        .value("water", molSys::atom_state_type::water)
        .value("interfacial", molSys::atom_state_type::interfacial)
        .value("clathrate", molSys::atom_state_type::clathrate)
        .value("interClathrate", molSys::atom_state_type::interClathrate)
        .value("unclassified", molSys::atom_state_type::unclassified)
        .value("reCubic", molSys::atom_state_type::reCubic)
        .value("reHex", molSys::atom_state_type::reHex);

    nb::class_<molSys::Result>(
        m, "Result", "Bond correlation result: classifier (bond type) and c_value.")
        .def(nb::init<>())
        .def_rw("classifier", &molSys::Result::classifier)
        .def_rw("c_value", &molSys::Result::c_value)
        .def("__repr__", [](const molSys::Result &self_C) {
            std::uintptr_t ptr_val = std::uintptr_t(&self_C);
            return std::format("<Result mem_loc:{:x}>", static_cast<uint64_t>(ptr_val));
        });

    nb::class_<molSys::PointCloud<molSys::Point<double>, double>>(
        m, "PointCloudDouble", "Collection of points for a single frame, with box dimensions.")
        .def(nb::init<>())
        .def_rw("pts", &molSys::PointCloud<molSys::Point<double>, double>::pts)
        .def_rw("currentFrame", &molSys::PointCloud<molSys::Point<double>, double>::currentFrame)
        .def_rw("nop", &molSys::PointCloud<molSys::Point<double>, double>::nop)
        .def_rw("box", &molSys::PointCloud<molSys::Point<double>, double>::box)
        .def_rw("boxLow", &molSys::PointCloud<molSys::Point<double>, double>::boxLow)
        .def_rw("idIndexMap", &molSys::PointCloud<molSys::Point<double>, double>::idIndexMap);

    nb::class_<chill::SteinhardtQl>(
        m,
        "SteinhardtQl",
        "Per-particle Steinhardt order parameters of a single degree l: "
        "the local ql and the neighbour-averaged qlBar.")
        .def(nb::init<>())
        .def_ro("ql", &chill::SteinhardtQl::ql)
        .def_ro("qlBar", &chill::SteinhardtQl::qlBar)
        .def("__repr__", [](const chill::SteinhardtQl &self_C) {
            std::uintptr_t ptr_val = std::uintptr_t(&self_C);
            return std::format("<SteinhardtQl mem_loc:{:x}>", static_cast<uint64_t>(ptr_val));
        });

    // I/O (lambdas hide the yCloud* in/out parameter, creating it internally)
    m.def(
        "readXYZ", &sinp::readXYZ, "Read atom coordinates from an XYZ file.", nb::arg("filename"));
    m.def(
        "readLammpsTrjreduced",
        [](std::string filename,
           int targetFrame,
           int typeI,
           bool isSlice,
           std::array<double, 3> coordLow,
           std::array<double, 3> coordHigh) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readLammpsTrjreduced(
                filename, targetFrame, yCloud, typeI, isSlice, coordLow, coordHigh);
        },
        "One LAMMPS frame, one type. isSlice drops atoms outside the box. "
        "nop is the kept count. An axis with lo == hi is unconstrained.",
        nb::arg("filename"),
        nb::arg("targetFrame"),
        nb::arg("typeI"),
        nb::arg("isSlice"),
        nb::arg("coordLow"),
        nb::arg("coordHigh"));
    m.def(
        "readLammpsTrjO",
        [](std::string filename,
           int targetFrame,
           int typeO,
           bool isSlice,
           std::array<double, 3> coordLow,
           std::array<double, 3> coordHigh) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readLammpsTrjO(
                filename, targetFrame, yCloud, typeO, isSlice, coordLow, coordHigh);
        },
        "One LAMMPS frame, one type (the O is historical). isSlice sets "
        "inSlice and does not drop. nop is the type-filtered count.",
        nb::arg("filename"),
        nb::arg("targetFrame"),
        nb::arg("typeO"),
        nb::arg("isSlice"),
        nb::arg("coordLow"),
        nb::arg("coordHigh"));
    m.def(
        "readLammpsTrj",
        [](std::string filename,
           int targetFrame,
           bool isSlice,
           std::array<double, 3> coordLow,
           std::array<double, 3> coordHigh) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readLammpsTrj(
                filename, targetFrame, yCloud, isSlice, coordLow, coordHigh);
        },
        "Read a LAMMPS trajectory frame with all atom types.",
        nb::arg("filename"),
        nb::arg("targetFrame"),
        nb::arg("isSlice"),
        nb::arg("coordLow"),
        nb::arg("coordHigh"));
    m.def("readBonds",
          &sinp::readBonds,
          "Read bond connectivity from a formatted bond file.",
          nb::arg("filename"));
    m.def("atomInSlice",
          &sinp::atomInSlice,
          "True when each component is in [lo, hi], or that axis has lo == hi.",
          nb::arg("x"),
          nb::arg("y"),
          nb::arg("z"),
          nb::arg("coordLow"),
          nb::arg("coordHigh"));
#ifdef SEAMS_HAS_CHEMFILES
    m.def(
        "readChemfiles",
        [](std::string filename, int targetFrame, int typeFilter) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readChemfiles(filename, targetFrame, yCloud, typeFilter);
        },
        "Read any trajectory format supported by chemfiles (PDB, GRO, DCD, etc.).",
        nb::arg("filename"),
        nb::arg("targetFrame"),
        nb::arg("typeFilter") = -1);
#endif
#ifdef SEAMS_HAS_READCON
    m.def(
        "readCon",
        [](std::string filename, int targetFrame) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readCon(filename, targetFrame, yCloud);
        },
        "Read a .con format file (eOn saddle point search trajectories).",
        nb::arg("filename"),
        nb::arg("targetFrame"));
#endif

    // Neighbours
    m.def("clearNeighbourList",
          &nneigh::clearNeighbourList,
          "Free memory for a neighbour list.",
          nb::arg("nList"));
    m.def("getNewNeighbourListByIndex",
          &nneigh::getNewNeighbourListByIndex,
          "Build a neighbour list by index using a distance cutoff.",
          nb::arg("yCloud"),
          nb::arg("cutoff"));
    m.def("mutualNearestUnlike",
          &nneigh::mutualNearestUnlike,
          "Mutual nearest unlike cloud-index pairs.",
          nb::arg("yCloud"),
          nb::arg("typeI"),
          nb::arg("typeJ"));
    nb::class_<clump::Domain>(
        m, "Domain", "Largest connected component statistics for a site mask.")
        .def_ro("largest", &clump::Domain::largest)
        .def_ro("subset", &clump::Domain::subset)
        .def_ro("percolation", &clump::Domain::percolation);
    m.def("largestDomain",
          &clump::largestDomain,
          "Largest connected component of a Boolean site mask.",
          nb::arg("yCloud"),
          nb::arg("nList"),
          nb::arg("mask"));
    using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
    using KnnTypeI = std::vector<std::vector<int>> (*)(
        const Cloud &, int, double, int, bool);
    using KnnTypes = std::vector<std::vector<int>> (*)(
        const Cloud &, int, double, const std::vector<int> &, bool);
    m.def("kNearestNeighbourList",
          static_cast<KnnTypeI>(&nneigh::kNearestNeighbourList),
          "Exact k-nearest bonded graph, union- or mutually-symmetrized "
          "(cell-list candidates with a brute-force fallback).",
          nb::arg("yCloud"),
          nb::arg("k"),
          nb::arg("candidateCutoff"),
          nb::arg("typeI"),
          nb::arg("mutual") = true);
    m.def("kNearestNeighbourList",
          static_cast<KnnTypes>(&nneigh::kNearestNeighbourList),
          "k-nearest bonded graph restricted to a type set.",
          nb::arg("yCloud"),
          nb::arg("k"),
          nb::arg("candidateCutoff"),
          nb::arg("types"),
          nb::arg("mutual") = true);
    m.def("shellSeparation",
          &nneigh::shellSeparation,
          "Certificate pair (max k-th distance, min (k+1)-th distance) for "
          "the exact reduction of the k-nearest graph to a cutoff graph.",
          nb::arg("yCloud"),
          nb::arg("k"),
          nb::arg("typeI"));
    m.def("halfNeighList",
          &nneigh::halfNeighList,
          "Build a half neighbour list (each pair stored once) for one atom type.",
          nb::arg("yCloud"),
          nb::arg("rcutoff"),
          nb::arg("typeI"));
    m.def("neighbourListByIndex",
          &nneigh::neighbourListByIndex,
          "Convert an atom-ID neighbour list to an index-based neighbour list.",
          nb::arg("yCloud"),
          nb::arg("nList"));
    m.def("neighList",
          &nneigh::neighList,
          "Build a full neighbour list for two atom types within a cutoff.",
          nb::arg("yCloud"),
          nb::arg("rcutoff"),
          nb::arg("typeI"),
          nb::arg("typeJ"));
    m.def("neighListO",
          &nneigh::neighListO,
          "Build a full neighbour list for a single atom type within a cutoff.",
          nb::arg("rcutoff"),
          nb::arg("yCloud"),
          nb::arg("typeI"));
    m.def("neighListPair",
          &nneigh::neighListPair,
          "I-J neighbour list. Like-type reuses neighListO; unlike-type "
          "pairs use the dump MIC.",
          nb::arg("rcutoff"),
          nb::arg("yCloud"),
          nb::arg("typeI"),
          nb::arg("typeJ"));

    // Bonds
    m.def("createBondsFromCages",
          &bond::createBondsFromCages,
          "Create bond connectivity from rings and cage information.",
          nb::arg("rings"),
          nb::arg("cageList"),
          nb::arg("type"),
          nb::arg("nRings"));
    m.def("getHbondDistanceOH",
          &bond::getHbondDistanceOH,
          "Compute the O-H hydrogen bond distance between two atoms.",
          nb::arg("oCloud"),
          nb::arg("hCloud"),
          nb::arg("oAtomIndex"),
          nb::arg("hAtomIndex"));
    m.def("populateHbonds",
          &bond::populateHbonds,
          "Build the hydrogen bond network from a trajectory and neighbour "
          "list. distCutoff/angleCutoff default to the water criterion "
          "(2.42 Angstrom, 30 degrees).",
          nb::arg("filename"),
          nb::arg("yCloud"),
          nb::arg("nList"),
          nb::arg("targetFrame"),
          nb::arg("Htype"),
          nb::arg("distCutoff") = 2.42,
          nb::arg("angleCutoff") = 30.0);
    m.def("populateHbondsWithInputClouds",
          &bond::populateHbondsWithInputClouds,
          "Build hydrogen bonds using pre-loaded oxygen and hydrogen point "
          "clouds. distCutoff/angleCutoff default to the water criterion "
          "(2.42 Angstrom, 30 degrees).",
          nb::arg("yCloud"),
          nb::arg("hCloud"),
          nb::arg("nList"),
          nb::arg("distCutoff") = 2.42,
          nb::arg("angleCutoff") = 30.0);
    m.def("donatedHydrogenBond",
          &bond::donatedHydrogenBond,
          "Geometric O-O-H test for one donor-acceptor assignment. "
          "donorHs are hCloud indices on the donor.",
          nb::arg("yCloud"),
          nb::arg("hCloud"),
          nb::arg("acceptorIndex"),
          nb::arg("donorIndex"),
          nb::arg("donorHs"),
          nb::arg("distCutoff") = 2.42,
          nb::arg("angleCutoff") = 30.0);
    m.def("populateHbondsFromDonors",
          &bond::populateHbondsFromDonors,
          "Hydrogen-bond network from an explicit donor-H set. donorHs "
          "are hCloud indices; each is paired with the heavy atom that "
          "shares its molID.",
          nb::arg("yCloud"),
          nb::arg("hCloud"),
          nb::arg("nList"),
          nb::arg("donorHs"),
          nb::arg("distCutoff") = 2.42,
          nb::arg("angleCutoff") = 30.0);
    m.def("trimBonds",
          &bond::trimBonds,
          "Remove duplicate bonds from a bond list.",
          nb::arg("bonds"));

    // Ring analysis (Franzblau)
    m.def("clearGraph",
          &primitive::clearGraph,
          "Free memory for a graph object.",
          nb::arg("currentGraph"));
    m.def("countAllRingsFromIndex",
          &primitive::countAllRingsFromIndex,
          "Find all possible rings (including non-shortest-path) up to maxDepth.",
          nb::arg("neighHbondList"),
          nb::arg("maxDepth"));
    m.def("ringNetwork",
          &primitive::ringNetwork,
          "Find all primitive (shortest-path) rings up to maxDepth.",
          nb::arg("nList"),
          nb::arg("maxDepth"));
    nb::class_<primitive::RingUpdater>(m, "RingUpdater")
        .def(nb::init<int>(), nb::arg("maxDepth"))
        .def("update",
             &primitive::RingUpdater::update,
             "Exact incremental primitive rings for this neighbour list.",
             nb::arg("nList"))
        .def("lastRecomputedSources", &primitive::RingUpdater::lastRecomputedSources)
        .def("lastBallsRefreshed", &primitive::RingUpdater::lastBallsRefreshed);
    m.def(
        "cageAffiliation",
        [](const std::vector<std::vector<int>> &rings,
           const std::vector<std::vector<int>> &nList) {
            const auto a = ring::cageAffiliation(rings, nList);
            return std::make_pair(a.hc, a.ddc);
        },
        "Order-free per-ring cage classification: (hc, ddc) flag vectors.",
        nb::arg("rings"),
        nb::arg("nList"));
    m.def(
        "findBySignature",
        [](const std::vector<std::vector<int>> &rings,
           const std::vector<std::vector<int>> &nList,
           const std::string &spec) {
            const auto sig = cage::Signature::parse(spec);
            const auto found = cage::findBySignature(rings, nList, sig);
            nb::list out;
            for (const auto &c : found) {
                nb::dict row;
                row["signature"] = c.signature.str();
                row["faces"] = c.faces;
                row["vertices"] = c.vertices;
                row["certificate"] = c.certificate;
                out.append(row);
            }
            return out;
        },
        "Closed polyhedra matching a ring-size census or named table entry.",
        nb::arg("rings"),
        nb::arg("nList"),
        nb::arg("spec"));
    nb::class_<ring::AffiliationUpdater>(m, "AffiliationUpdater")
        .def(nb::init<>())
        .def(
            "update",
            [](ring::AffiliationUpdater &self,
               const std::vector<std::vector<int>> &rings,
               const std::vector<std::vector<int>> &nList) {
                const auto &a = self.update(rings, nList);
                return std::make_pair(a.hc, a.ddc);
            },
            "Exact incremental per-ring cage classification for this frame.",
            nb::arg("rings"),
            nb::arg("nList"))
        .def("lastReclassified", &ring::AffiliationUpdater::lastReclassified);
    m.def(
        "seededCageAffiliation",
        [](const std::vector<std::vector<int>> &strictRings,
           const std::vector<std::vector<int>> &strictNList,
           const std::vector<std::vector<int>> &permissiveRings,
           const std::vector<std::vector<int>> &permissiveNList,
           bool ringAdjacentCompletion) {
            const auto a = ring::seededCageAffiliation(strictRings,
                                                       strictNList,
                                                       permissiveRings,
                                                       permissiveNList,
                                                       ringAdjacentCompletion);
            return std::make_pair(a.hc, a.ddc);
        },
        "Seeded (hysteresis) per-atom cage flags: strict-graph seeds, "
        "permissive-graph completion, component-gated acceptance. With "
        "ringAdjacentCompletion a permissive six-ring whose vertices all carry a label "
        "but one fills that last vertex, repeated to a fixed point.",
        nb::arg("strictRings"),
        nb::arg("strictNList"),
        nb::arg("permissiveRings"),
        nb::arg("permissiveNList"),
        nb::arg("ringAdjacentCompletion") = false);
    m.def("populateGraphFromIndices",
          &primitive::populateGraphFromIndices,
          "Create a graph object from an index-based neighbour list.",
          nb::arg("nList"));
    m.def("populateGraphFromNListID",
          &primitive::populateGraphFromNListID,
          "Create a graph object from an atom-ID neighbour list and point cloud.",
          nb::arg("yCloud"),
          nb::arg("neighHbondList"));
    m.def("removeNonSPrings",
          &primitive::removeNonSPrings,
          "Remove non-shortest-path rings using the Franzblau criterion.",
          nb::arg("fullGraph"));
    m.def("restoreEdgesFromIndices",
          &primitive::restoreEdgesFromIndices,
          "Restore graph edges from an index-based neighbour list.",
          nb::arg("fullGraph"),
          nb::arg("nList"));

    // Ring classification
    m.def("assignPolygonType",
          &ring::assignPolygonType,
          "Assign atom types based on the ring size of n-membered rings.",
          nb::arg("rings"),
          nb::arg("atomTypes"),
          nb::arg("nRings"));
    m.def("assignPrismType",
          &ring::assignPrismType,
          "Assign atom types for atoms belonging to prism rings.",
          nb::arg("rings"),
          nb::arg("listPrism"),
          nb::arg("ringSize"),
          nb::arg("ringType"),
          nb::arg("atomTypes"),
          nb::arg("atomState"));
    m.def("clearRingList",
          &ring::clearRingList,
          "Free memory for a list of rings.",
          nb::arg("rings"));
    m.def("compareRings",
          &ring::compareRings,
          "Check whether two unordered rings contain the same elements.",
          nb::arg("ring1"),
          nb::arg("ring2"));
    m.def("commonElementsInThreeRings",
          &ring::commonElementsInThreeRings,
          "Check whether three rings share at least one common element.",
          nb::arg("ring1"),
          nb::arg("ring2"),
          nb::arg("ring3"));
    m.def("deformedPrismTypes",
          &ring::deformedPrismTypes,
          "Get atom type values for deformed prisms.",
          nb::arg("atomState"),
          nb::arg("atomTypes"),
          nb::arg("maxDepth"));
    m.def("discardExtraTetragonBlocks",
          &ring::discardExtraTetragonBlocks,
          "Discard duplicate 4-membered ring pairs that are parallel in one dimension.",
          nb::arg("basal1"),
          nb::arg("basal2"),
          nb::arg("yCloud"));
    m.def("findPrisms",
          &ring::findPrisms,
          "Identify which rings form prism blocks.",
          nb::arg("rings"),
          nb::arg("ringType"),
          nb::arg("nPerfectPrisms"),
          nb::arg("nImperfectPrisms"),
          nb::arg("nList"),
          nb::arg("rmsdPerAtom"),
          nb::arg("doShapeMatching"),
          nb::arg("yCloud"));
    m.def("findsCommonElements",
          &ring::findsCommonElements,
          "Return the common elements shared by two rings.",
          nb::arg("ring1"),
          nb::arg("ring2"));
    m.def("findTripletInRing",
          &ring::findTripletInRing,
          "Search for a triplet of atoms within a ring.",
          nb::arg("ring"),
          nb::arg("triplet"));
    m.def("getSingleRingSize",
          &ring::getSingleRingSize,
          "Extract rings of a specific size from a list of all rings.",
          nb::arg("rings"),
          nb::arg("ringSize"));
    m.def("hasCommonElements",
          &ring::hasCommonElements,
          "Check whether two rings share any common elements.",
          nb::arg("ring1"),
          nb::arg("ring2"));
    m.def("basalPrismConditions",
          &ring::basalPrismConditions,
          "Test whether two rings satisfy strict basal prism conditions.",
          nb::arg("nList"),
          nb::arg("basal1"),
          nb::arg("basal2"));
    m.def("relaxedPrismConditions",
          &ring::relaxedPrismConditions,
          "Test whether two rings satisfy relaxed prism conditions (at least one bond).",
          nb::arg("nList"),
          nb::arg("basal1"),
          nb::arg("basal2"));

    // Topology analysis
    m.def("polygonRingAnalysis",
          &ring::polygonRingAnalysis,
          "Classify rings in a quasi-2D monolayer system and write output.",
          nb::arg("path"),
          nb::arg("rings"),
          nb::arg("nList"),
          nb::arg("yCloud"),
          nb::arg("maxDepth"),
          nb::arg("sheetArea"),
          nb::arg("firstFrame"));
    m.def("bulkPolygonRingAnalysis",
          &ring::bulkPolygonRingAnalysis,
          "Classify rings in a bulk system and write output.",
          nb::arg("path"),
          nb::arg("rings"),
          nb::arg("nList"),
          nb::arg("yCloud"),
          nb::arg("maxDepth"),
          nb::arg("firstFrame"));
    m.def("prismAnalysis",
          &ring::prismAnalysis,
          "Run prism identification on rings up to maxDepth and write output.",
          nb::arg("path"),
          nb::arg("rings"),
          nb::arg("nList"),
          nb::arg("yCloud"),
          nb::arg("maxDepth"),
          nb::arg("atomID"),
          nb::arg("firstFrame"),
          nb::arg("currentFrame"),
          nb::arg("doShapeMatching"));
    m.def("rmAxialTranslations",
          &ring::rmAxialTranslations,
          "Remove axial translations from an ice nanotube for visualization.",
          nb::arg("yCloud"),
          nb::arg("atomID"),
          nb::arg("firstFrame"),
          nb::arg("currentFrame"));

    // TUM (Topological Unit Matching)
    m.def("atomsFromCages",
          &tum3::atomsFromCages,
          "Get atom indices belonging to cages in a given cluster.",
          nb::arg("rings"),
          nb::arg("cageList"),
          nb::arg("clusterCages"));
    m.def("averageRMSDatom",
          &tum3::averageRMSDatom,
          "Average the per-atom RMSD over the number of shared cages.",
          nb::arg("rmsdPerAtom"),
          nb::arg("noOfCommonAtoms"));
    m.def("buildRefDDC",
          &tum3::buildRefDDC,
          "Build a reference double-diamond cage from a template XYZ file.",
          nb::arg("fileName"));
    m.def("buildRefHC",
          &tum3::buildRefHC,
          "Build a reference hexagonal cage from a template XYZ file.",
          nb::arg("fileName"));
    m.def("clusterCages",
          &tum3::clusterCages,
          "Cluster cages using Stillinger's algorithm and write XYZ output.",
          nb::arg("yCloud"),
          nb::arg("path"),
          nb::arg("rings"),
          nb::arg("cageList"),
          nb::arg("numHC"),
          nb::arg("numDDC"));
    m.def("shapeMatchDDC",
          &tum3::shapeMatchDDC,
          "Shape-match a target double-diamond cage against a reference.",
          nb::arg("yCloud"),
          nb::arg("refPoints"),
          nb::arg("cageList"),
          nb::arg("cageIndex"),
          nb::arg("rings"),
          nb::arg("quat"),
          nb::arg("rmsd"));
    m.def("shapeMatchHC",
          &tum3::shapeMatchHC,
          "Shape-match a target hexagonal cage against a reference.",
          nb::arg("yCloud"),
          nb::arg("refPoints"),
          nb::arg("cageUnit"),
          nb::arg("rings"),
          nb::arg("nList"),
          nb::arg("quat"),
          nb::arg("rmsd"));
    m.def("topoBulkCriteria",
          &tum3::topoBulkCriteria,
          "Find HCs and DDCs in a bulk system using topological criteria.",
          nb::arg("path"),
          nb::arg("rings"),
          nb::arg("nList"),
          nb::arg("yCloud"),
          nb::arg("firstFrame"),
          nb::arg("numHC"),
          nb::arg("numDDC"),
          nb::arg("ringType"));
    m.def("topoUnitMatchingBulk",
          &tum3::topoUnitMatchingBulk,
          "Run full topological unit matching for bulk water.",
          nb::arg("path"),
          nb::arg("rings"),
          nb::arg("nList"),
          nb::arg("yCloud"),
          nb::arg("firstFrame"),
          nb::arg("printClusters"),
          nb::arg("onlyTetrahedral"),
          nb::arg("templatePath") = "templates");
    m.def("updateRMSDatom",
          &tum3::updateRMSDatom,
          "Update per-atom RMSD from a cage shape-matching result.",
          nb::arg("rings"),
          nb::arg("cageUnit"),
          nb::arg("rmsd"),
          nb::arg("rmsdPerAtom"),
          nb::arg("noOfCommonAtoms"),
          nb::arg("atomTypes"));

    // Selection
    m.def("getPointCloudOneAtomType",
          &gen::getPointCloudOneAtomType,
          "Extract a point cloud containing only atoms of a given type.",
          nb::arg("yCloud"),
          nb::arg("outCloud"),
          nb::arg("atomTypeI"),
          nb::arg("isSlice"),
          nb::arg("coordLow"),
          nb::arg("coordHigh"));
    m.def("atomsInSingleSlice",
          &gen::atomsInSingleSlice,
          "Mark atoms inside a rectangular volume slice.",
          nb::arg("yCloud"),
          nb::arg("clearPreviousSliceSelection"),
          nb::arg("coordLow"),
          nb::arg("coordHigh"));
    m.def("moleculesInSingleSlice",
          &gen::moleculesInSingleSlice,
          "Mark whole molecules as in-slice if any atom is inside the region.",
          nb::arg("yCloud"),
          nb::arg("clearPreviousSliceSelection"),
          nb::arg("coordLow"),
          nb::arg("coordHigh"));
    m.def("setAtomsWithSameMolID",
          &gen::setAtomsWithSameMolID,
          "Set the inSlice flag for all atoms sharing a given molecule ID.",
          nb::arg("yCloud"),
          nb::arg("molIDAtomIDmap"),
          nb::arg("molID"),
          nb::arg("inSliceValue"));
    m.def("getEdgeMoleculesInRings",
          &ring::getEdgeMoleculesInRings,
          "Select edge molecules in rings that straddle the slice boundary.",
          nb::arg("rings"),
          nb::arg("oCloud"),
          nb::arg("yCloud"),
          nb::arg("identicalCloud"),
          nb::arg("coordLow"),
          nb::arg("coordHigh"));
    m.def("printSliceGetEdgeMoleculesInRings",
          &ring::printSliceGetEdgeMoleculesInRings,
          "Select edge molecules in rings and write slice output files.",
          nb::arg("path"),
          nb::arg("rings"),
          nb::arg("oCloud"),
          nb::arg("yCloud"),
          nb::arg("coordLow"),
          nb::arg("coordHigh"),
          nb::arg("identicalCloud"));

    // Bond-classification rule sets: CHILL and CHILL+ are two registered
    // instances of the same engine, and other materials register their own
    nb::class_<chill::BondClassifier>(m, "BondClassifier")
        .def(nb::init<>())
        .def(
            "__init__",
            [](chill::BondClassifier *self,
               double staggeredMax,
               double eclipsedMin,
               double eclipsedMax,
               int coordinationNumber) {
                new (self) chill::BondClassifier{
                    staggeredMax, eclipsedMin, eclipsedMax, coordinationNumber};
            },
            nb::arg("staggeredMax"),
            nb::arg("eclipsedMin"),
            nb::arg("eclipsedMax"),
            nb::arg("coordinationNumber") = 4)
        .def_rw("staggeredMax", &chill::BondClassifier::staggeredMax)
        .def_rw("eclipsedMin", &chill::BondClassifier::eclipsedMin)
        .def_rw("eclipsedMax", &chill::BondClassifier::eclipsedMax)
        .def_rw("coordinationNumber", &chill::BondClassifier::coordinationNumber);
    m.def("chillRule", &chill::chillRule, "The CHILL water rule set.");
    m.def("chillPlusRule", &chill::chillPlusRule, "The CHILL+ water rule set.");
    m.def("bondClassifier",
          &chill::bondClassifier,
          "Look up a registered bond-classification rule set by name.",
          nb::arg("name"));
    m.def("registerBondClassifier",
          &chill::registerBondClassifier,
          "Register (or replace) a named bond-classification rule set.",
          nb::arg("name"),
          nb::arg("rule"));
    m.def("bondClassifierNames",
          &chill::bondClassifierNames,
          "Names of every registered bond-classification rule set.");
    m.def(
        "classifyBonds",
        [](molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           const std::vector<std::vector<int>> &nList,
           const chill::BondClassifier &rule,
           bool isSlice) -> molSys::PointCloud<molSys::Point<double>, double> & {
            chill::classifyBonds(yCloud, nList, rule, isSlice);
            return yCloud;
        },
        "Compute and classify bond correlations under an arbitrary rule set.",
        nb::arg("yCloud"),
        nb::arg("nList"),
        nb::arg("rule"),
        nb::arg("isSlice") = false,
        nb::rv_policy::reference);

    // CHILL/CHILL+ classification. The C++ routines mutate yCloud in place
    // and return void. The Python bindings return the same object so that
    // `cloud = yoda.getCorrel(yCloud=cloud, ...)` keeps the caller's
    // reference, which is the contract the 2.0 bindings shipped.
    m.def(
        "getCorrelPlus",
        [](molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           const std::vector<std::vector<int>> &nList,
           bool isSlice,
           int coordinationNumber) -> molSys::PointCloud<molSys::Point<double>, double> & {
            chill::getCorrelPlus(yCloud, nList, isSlice, coordinationNumber);
            return yCloud;
        },
        "Compute CHILL+ bond-order correlations and classify bond types. "
        "coordinationNumber=4 is the validated water scheme; non-positive "
        "uses each atom's whole neighbour row.",
        nb::arg("yCloud"),
        nb::arg("nList"),
        nb::arg("isSlice"),
        nb::arg("coordinationNumber") = 4,
        nb::rv_policy::reference);
    m.def(
        "getIceTypePlus",
        [](molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           const std::vector<std::vector<int>> &nList,
           std::string path,
           int firstFrame,
           bool isSlice,
           std::string outputFileName) -> molSys::PointCloud<molSys::Point<double>, double> & {
            chill::getIceTypePlus(yCloud, nList, path, firstFrame, isSlice, outputFileName);
            return yCloud;
        },
        "Classify each atom's ice type using CHILL+ and write to file.",
        nb::arg("yCloud"),
        nb::arg("nList"),
        nb::arg("path"),
        nb::arg("firstFrame"),
        nb::arg("isSlice"),
        nb::arg("outputFileName"),
        nb::rv_policy::reference);
    m.def(
        "getIceTypePlusNoPrint",
        [](molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           const std::vector<std::vector<int>> &nList,
           bool isSlice) -> molSys::PointCloud<molSys::Point<double>, double> & {
            chill::getIceTypePlusNoPrint(yCloud, nList, isSlice);
            return yCloud;
        },
        "Classify each atom's ice type using CHILL+. Does not write a file.",
        nb::arg("yCloud"),
        nb::arg("nList"),
        nb::arg("isSlice") = false,
        nb::rv_policy::reference);
    m.def(
        "getIceTypeNoPrint",
        [](molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           const std::vector<std::vector<int>> &nList,
           bool isSlice) -> molSys::PointCloud<molSys::Point<double>, double> & {
            chill::getIceTypeNoPrint(yCloud, nList, isSlice);
            return yCloud;
        },
        "Classify each atom's ice type using CHILL. Does not write a file.",
        nb::arg("yCloud"),
        nb::arg("nList"),
        nb::arg("isSlice") = false,
        nb::rv_policy::reference);
    m.def(
        "getCorrel",
        [](molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           const std::vector<std::vector<int>> &nList,
           bool isSlice,
           int coordinationNumber) -> molSys::PointCloud<molSys::Point<double>, double> & {
            chill::getCorrel(yCloud, nList, isSlice, coordinationNumber);
            return yCloud;
        },
        "Compute CHILL bond-order correlations and classify bond types. "
        "coordinationNumber=4 is the validated water scheme; non-positive "
        "uses each atom's whole neighbour row.",
        nb::arg("yCloud"),
        nb::arg("nList"),
        nb::arg("isSlice"),
        nb::arg("coordinationNumber") = 4,
        nb::rv_policy::reference);
    m.def(
        "getIceType",
        [](molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           const std::vector<std::vector<int>> &nList,
           std::string path,
           int firstFrame,
           bool isSlice,
           std::string outputFileName) -> molSys::PointCloud<molSys::Point<double>, double> & {
            chill::getIceType(yCloud, nList, path, firstFrame, isSlice, outputFileName);
            return yCloud;
        },
        "Classify each atom's ice type using CHILL and write to file.",
        nb::arg("yCloud"),
        nb::arg("nList"),
        nb::arg("path"),
        nb::arg("firstFrame"),
        nb::arg("isSlice"),
        nb::arg("outputFileName"),
        nb::rv_policy::reference);
    m.def("getq6",
          &chill::getq6,
          "Compute the q6 bond order parameter for all atoms.",
          nb::arg("yCloud"),
          nb::arg("nList"),
          nb::arg("isSlice"));
    m.def(
        "reclassifyWater",
        [](molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           std::vector<double> &q6) -> molSys::PointCloud<molSys::Point<double>, double> & {
            chill::reclassifyWater(yCloud, q6);
            return yCloud;
        },
        "Reclassify water molecules using averaged q6 and q3 parameters.",
        nb::arg("yCloud"),
        nb::arg("q6"),
        nb::rv_policy::reference);
    m.def("printIceType",
          &chill::printIceType,
          "Print the ice type classification for the current frame.",
          nb::arg("yCloud"),
          nb::arg("path"),
          nb::arg("firstFrame"),
          nb::arg("isSlice"),
          nb::arg("outputFileName"));
    m.def("steinhardtQl",
          &chill::steinhardtQl,
          "Compute the local and neighbour-averaged Steinhardt parameters "
          "of degree orderL (3, 4 or 6) for every particle.",
          nb::arg("yCloud"),
          nb::arg("nList"),
          nb::arg("orderL"));
    nb::enum_<chill::CrystalKind>(m, "CrystalKind")
        .value("other", chill::CrystalKind::other)
        .value("sc", chill::CrystalKind::sc)
        .value("fcc", chill::CrystalKind::fcc)
        .value("hcp", chill::CrystalKind::hcp)
        .value("bcc", chill::CrystalKind::bcc);
    nb::class_<chill::TemplateHit>(m, "TemplateHit")
        .def_ro("kind", &chill::TemplateHit::kind)
        .def_ro("rmsd", &chill::TemplateHit::rmsd)
        .def_prop_ro("name", [](const chill::TemplateHit &h) { return std::string(h.name); });
    nb::class_<chill::LinearClassifier>(m, "LinearClassifier")
        .def(nb::init<>())
        .def_rw("ridge", &chill::LinearClassifier::ridge)
        .def_rw("labels", &chill::LinearClassifier::labels)
        .def_ro("n_classes", &chill::LinearClassifier::nClasses)
        .def_ro("n_feat", &chill::LinearClassifier::nFeat)
        .def("fit", &chill::LinearClassifier::fit, nb::arg("X"), nb::arg("y"))
        .def("predict", &chill::LinearClassifier::predict, nb::arg("x"));
    m.def("steinhardtQlVoronoi",
          &chill::steinhardtQlVoronoi,
          "Voronoi facet-area weighted Steinhardt parameters.",
          nb::arg("yCloud"),
          nb::arg("candidateCutoff"),
          nb::arg("orderL"));
    m.def("classifyTemplates",
          &chill::classifyTemplates,
          "IRA/Horn overlay onto FCC, HCP, BCC and SC neighbour shells.",
          nb::arg("yCloud"),
          nb::arg("nList"),
          nb::arg("kNeigh") = 12);
    m.def("soapSpectrum",
          &chill::soapSpectrum,
          "SOAP power spectrum of one particle.",
          nb::arg("yCloud"),
          nb::arg("iatom"),
          nb::arg("nList"),
          nb::arg("nMax"),
          nb::arg("lMax"),
          nb::arg("rcut"));
    m.def("soapSpectrumAll",
          &chill::soapSpectrumAll,
          "SOAP power spectrum of every particle.",
          nb::arg("yCloud"),
          nb::arg("nList"),
          nb::arg("nMax"),
          nb::arg("lMax"),
          nb::arg("rcut"));
    m.def("voronoiFeature",
          &chill::voronoiFeature,
          "Per-atom [q4, q6, q8] from the Voronoi-weighted Steinhardt path.",
          nb::arg("yCloud"),
          nb::arg("iatom"),
          nb::arg("candidateCutoff"));
    m.def("voronoiFeatures",
          &chill::voronoiFeatures,
          "[q4, q6, q8] for every particle from one Voronoi pass per order.",
          nb::arg("yCloud"),
          nb::arg("candidateCutoff"));

    m.def("ira_available", &ira::available, "True when this build linked libira (IRA/SOFI).");
#ifdef SEAMS_HAS_IRA
    m.def(
        "ira_match",
        [](const std::vector<std::vector<double>> &refPts,
           const std::vector<std::vector<double>> &tgtPts) {
            auto toMat = [](const std::vector<std::vector<double>> &rows) {
                Eigen::MatrixXd m(static_cast<int>(rows.size()), 3);
                for (int i = 0; i < static_cast<int>(rows.size()); i++) {
                    m(i, 0) = rows[static_cast<size_t>(i)][0];
                    m(i, 1) = rows[static_cast<size_t>(i)][1];
                    m(i, 2) = rows[static_cast<size_t>(i)][2];
                }
                return m;
            };
            ira::Match out;
            const int err = ira::match(toMat(refPts), toMat(tgtPts), out);
            return nb::make_tuple(err, out.rmsd, out.hausdorff, out.quat, out.assignment);
        },
        "IRA overlay of two n x 3 point sets. Returns err, rmsd, hausdorff, "
        "quat, assignment.",
        nb::arg("ref"),
        nb::arg("target"));
    m.def(
        "sofi_point_group",
        [](const std::vector<std::vector<double>> &pts) {
            Eigen::MatrixXd m(static_cast<int>(pts.size()), 3);
            for (int i = 0; i < static_cast<int>(pts.size()); i++) {
                m(i, 0) = pts[static_cast<size_t>(i)][0];
                m(i, 1) = pts[static_cast<size_t>(i)][1];
                m(i, 2) = pts[static_cast<size_t>(i)][2];
            }
            ira::PointGroup pg;
            const int err = ira::pointGroup(m, pg);
            return nb::make_tuple(err, pg.symbol, pg.nOperations);
        },
        "SOFI point group of an n x 3 cloud. Returns err, symbol, n_ops.",
        nb::arg("coords"));
#endif

    // Spherical harmonics lookup tables
    m.def("lookupTableQ4Vec",
          &sph::lookupTableQ4Vec,
          "Lookup table for Q4 (m=0 to m=8).",
          nb::arg("angles"));
    m.def("lookupTableQ4",
          &sph::lookupTableQ4,
          "Lookup table for Q4 at a single m (m=0 to m=8).",
          nb::arg("m"),
          nb::arg("angles"));
    m.def("lookupTableQ8Vec",
          &sph::lookupTableQ8Vec,
          "Lookup table for Q8 (m=0 to m=16).",
          nb::arg("angles"));
    m.def("lookupTableQ8",
          &sph::lookupTableQ8,
          "Lookup table for Q8 at a single m (m=0 to m=16).",
          nb::arg("m"),
          nb::arg("angles"));

    // Voronoi facet-area weighted parameters
    nb::class_<chill::VoronoiWeights>(
        m,
        "VoronoiWeights",
        "Facet-sharing neighbours of one particle and their facet-area "
        "weights, normalised to sum to one.")
        .def(nb::init<>())
        .def_ro("neighbours", &chill::VoronoiWeights::neighbours)
        .def_ro("weights", &chill::VoronoiWeights::weights)
        .def_ro("certified", &chill::VoronoiWeights::certified);
    m.def("voronoiFacetWeights",
          &chill::voronoiFacetWeights,
          "Voronoi facet neighbours and area weights for every particle.",
          nb::arg("yCloud"),
          nb::arg("candidateCutoff"));

    // Output
    m.def("writeDump",
          &sout::writeDump,
          "Write a LAMMPS dump file for the current point cloud.",
          nb::arg("yCloud"),
          nb::arg("path"),
          nb::arg("outFile"));

    // Clustering
    m.def("clusterAnalysis",
          &clump::clusterAnalysis,
          "Cluster ice-like particles and return the largest ice cluster.",
          nb::arg("path"),
          nb::arg("iceCloud"),
          nb::arg("yCloud"),
          nb::arg("nList"),
          nb::arg("iceNeighbourList"),
          nb::arg("cutoff"),
          nb::arg("firstFrame"),
          nb::arg("bopAnalysis"));
    m.def("recenterClusterCloud",
          &clump::recenterClusterCloud,
          "Recenter a cluster point cloud for visualization.",
          nb::arg("iceCloud"),
          nb::arg("nList"));

    // RDF
    m.def("rdf2Danalysis_AA",
          &rdf2::rdf2Danalysis_AA,
          "Compute the 2D radial distribution function for identical atom types.",
          nb::arg("path"),
          nb::arg("rdfValues"),
          nb::arg("yCloud"),
          nb::arg("cutoff"),
          nb::arg("binwidth"),
          nb::arg("firstFrame"),
          nb::arg("finalFrame"));
    m.def(
        "partialRdf",
        [](const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           int typeI,
           int typeJ,
           double rmax,
           int nbins) {
            const auto h = rdf::partialRdf(yCloud, typeI, typeJ, rmax, nbins);
            return std::pair<std::vector<double>, std::vector<double>>(h.r, h.g);
        },
        "Partial 3D radial distribution function g_IJ(r). Returns (r, g).",
        nb::arg("yCloud"),
        nb::arg("typeI"),
        nb::arg("typeJ"),
        nb::arg("rmax"),
        nb::arg("nbins"));
    nb::class_<rdf::PartialRdf>(m,
                                "PartialRdf",
                                "Histogram from rdf::partialRdf: bin centres, g_IJ, pair counts, "
                                "and the dump-cell volume used to normalize.")
        .def_ro("r", &rdf::PartialRdf::r)
        .def_ro("g", &rdf::PartialRdf::g)
        .def_ro("count", &rdf::PartialRdf::count)
        .def_ro("rmax", &rdf::PartialRdf::rmax)
        .def_ro("binwidth", &rdf::PartialRdf::binwidth)
        .def_ro("volume", &rdf::PartialRdf::volume)
        .def_ro("typeI", &rdf::PartialRdf::typeI)
        .def_ro("typeJ", &rdf::PartialRdf::typeJ)
        .def_ro("nI", &rdf::PartialRdf::nI)
        .def_ro("nJ", &rdf::PartialRdf::nJ);
    m.def("partialRdfHist",
          &rdf::partialRdf,
          "Partial 3D RDF as a PartialRdf (r, g, count, volume, nI, nJ).",
          nb::arg("yCloud"),
          nb::arg("typeI"),
          nb::arg("typeJ"),
          nb::arg("rmax"),
          nb::arg("nbins"));
    m.def(
        "runningCN",
        [](const rdf::PartialRdf &h, std::optional<double> rhoJ) {
            const double rho
                = rhoJ.value_or((h.volume > 0.0) ? static_cast<double>(h.nJ) / h.volume : 0.0);
            return rdf::runningCN(h, rho);
        },
        "Running site-site CN: 4 pi rho_J int s^2 g(s) ds. "
        "rhoJ defaults to nJ / volume.",
        nb::arg("h"),
        nb::arg("rhoJ") = nb::none());
    m.def("firstMinimumBin",
          &rdf::firstMinimumBin,
          "Bin index of the first minimum of g after the first peak, or -1.",
          nb::arg("h"));
    m.def(
        "coordinationNumber",
        [](const rdf::PartialRdf &h, double rMax, std::optional<double> rhoJ) {
            const double rho
                = rhoJ.value_or((h.volume > 0.0) ? static_cast<double>(h.nJ) / h.volume : 0.0);
            return rdf::coordinationNumber(h, rMax, rho);
        },
        "Site-site CN integrated to rMax. rhoJ defaults to nJ / volume.",
        nb::arg("h"),
        nb::arg("rMax"),
        nb::arg("rhoJ") = nb::none());

    nb::enum_<site::Kind>(m, "SiteKind", "Per-atom site chemistry. Type 1 is not a kind.")
        .value("unspecified", site::Kind::unspecified)
        .value("cationHead", site::Kind::cationHead)
        .value("anion", site::Kind::anion)
        .value("tail", site::Kind::tail)
        .value("donorH", site::Kind::donorH)
        .value("acceptor", site::Kind::acceptor)
        .value("polar", site::Kind::polar)
        .value("apolar", site::Kind::apolar)
        .value("waterO", site::Kind::waterO)
        .value("waterH", site::Kind::waterH)
        .value("solvent", site::Kind::solvent);
    nb::enum_<site::Family>(m, "SiteFamily", "Input family; not inferred from types.")
        .value("waterIce", site::Family::waterIce)
        .value("ionicLiquid", site::Family::ionicLiquid)
        .value("moltenSalt", site::Family::moltenSalt)
        .value("des", site::Family::des)
        .value("electrolyte", site::Family::electrolyte)
        .value("confinedIL", site::Family::confinedIL)
        .value("confinedWater", site::Family::confinedWater)
        .value("networkFormer", site::Family::networkFormer);
    nb::class_<site::Table>(
        m, "SiteTable", "Map LAMMPS types (and optional atom-ID overrides) onto Kind.")
        .def(nb::init<>())
        .def_rw("family", &site::Table::family)
        .def_rw("typeToKind", &site::Table::typeToKind)
        .def_rw("atomOverride", &site::Table::atomOverride)
        .def("of", &site::Table::of, nb::arg("p"))
        .def("ofType", &site::Table::ofType, nb::arg("typeId"));
    nb::enum_<site::IonState>(m, "IonState", "An ion by its first water shell.")
        .value("liquid", site::IonState::liquid)
        .value("front", site::IonState::front)
        .value("ice", site::IonState::ice);
    nb::class_<site::IonEnvironment>(
        m, "IonEnvironment", "Ions read against a per-atom ice assignment.")
        .def_ro("ion", &site::IonEnvironment::ion)
        .def_ro("shell", &site::IonEnvironment::shell)
        .def_ro("iceFraction", &site::IonEnvironment::iceFraction)
        .def_ro("state", &site::IonEnvironment::state)
        .def_ro("members", &site::IonEnvironment::members)
        .def_ro("nIce", &site::IonEnvironment::nIce)
        .def_ro("nFront", &site::IonEnvironment::nFront)
        .def_ro("nLiquid", &site::IonEnvironment::nLiquid);
    m.def("ionEnvironment",
          &site::ionEnvironment,
          "Class each ion by its first water shell against iceFlag: every "
          "shell molecule labelled is ice, none is liquid, otherwise front.",
          nb::arg("yCloud"),
          nb::arg("iceFlag"),
          nb::arg("ionIndices"),
          nb::arg("waterType") = 1,
          nb::arg("cutoff") = 3.5);
    m.def("shellRingCensus",
          &site::shellRingCensus,
          "Census by size of the rings with at least one vertex in shell, up to "
          "maxRingSize: the rings of the water network that pass through an ion's "
          "first shell.",
          nb::arg("rings"),
          nb::arg("shell"),
          nb::arg("maxRingSize") = 7);
    nb::class_<topo::LocalKey>(m, "LocalKey", "Key of one atom's rooted bonded neighbourhood.")
        .def_ro("key", &topo::LocalKey::key)
        .def_ro("method", &topo::LocalKey::method)
        .def_ro("vertices", &topo::LocalKey::vertices)
        .def_ro("edges", &topo::LocalKey::edges);
    nb::class_<topo::FrameFingerprint>(
        m,
        "FrameFingerprint",
        "Local keys of every atom, their histogram, the ring census, one frame key.")
        .def_ro("key", &topo::FrameFingerprint::key)
        .def_ro("method", &topo::FrameFingerprint::method)
        .def_ro("atomKeys", &topo::FrameFingerprint::atomKeys)
        .def_ro("classes", &topo::FrameFingerprint::classes)
        .def_ro("ringCensus", &topo::FrameFingerprint::ringCensus)
        .def_ro("hops", &topo::FrameFingerprint::hops)
        .def_ro("coloured", &topo::FrameFingerprint::coloured);
    nb::class_<topo::KeyLibrary>(
        m, "KeyLibrary", "Local keys of reference structures with their labels.")
        .def(nb::init<>())
        .def_ro("method", &topo::KeyLibrary::method)
        .def_ro("hops", &topo::KeyLibrary::hops)
        .def_ro("coloured", &topo::KeyLibrary::coloured)
        .def_ro("labelOf", &topo::KeyLibrary::labelOf);
    nb::class_<topo::LibraryMatch>(m, "LibraryMatch", "Atoms named by a key library.")
        .def_ro("labels", &topo::LibraryMatch::labels)
        .def_ro("counts", &topo::LibraryMatch::counts)
        .def_ro("depth", &topo::LibraryMatch::depth)
        .def_ro("matched", &topo::LibraryMatch::matched);
    m.def("addToLibrary",
          &topo::addToLibrary,
          "Add every distinct key of a fingerprint under a label.",
          nb::arg("library"),
          nb::arg("fingerprint"),
          nb::arg("label"));
    m.def("writeLibrary", &topo::writeLibrary, "Text form of a key library.", nb::arg("library"));
    m.def("readLibrary", &topo::readLibrary, "Key library from its text form.", nb::arg("text"));
    m.def("matchLibrary",
          &topo::matchLibrary,
          "Name every atom of a fingerprint by a key library (same method, hops and "
          "colouring).",
          nb::arg("fingerprint"),
          nb::arg("library"));
    m.def("matchLibraries",
          &topo::matchLibraries,
          "Name every atom by key libraries at several hop counts: the deepest library "
          "that holds an atom's key names it, and depth records which one did.",
          nb::arg("rows"),
          nb::arg("libraries"),
          nb::arg("maxRingSize") = 7,
          nb::arg("colours") = std::vector<int>{});
    nb::class_<site::GuestOccupancy>(
        m, "GuestOccupancy", "Guests placed in enumerated cages by their periodic centroids.")
        .def_ro("guestsPerCage", &site::GuestOccupancy::guestsPerCage)
        .def_ro("cageOfGuest", &site::GuestOccupancy::cageOfGuest)
        .def_ro("centreDistance", &site::GuestOccupancy::centreDistance)
        .def_ro("occupied", &site::GuestOccupancy::occupied)
        .def_ro("multiply", &site::GuestOccupancy::multiply)
        .def_ro("free", &site::GuestOccupancy::free);
    m.def("guestOccupancy",
          &site::guestOccupancy,
          "Assign each guest to the nearest cage centroid within radius; cages are "
          "vertex index lists into the cloud.",
          nb::arg("yCloud"),
          nb::arg("cages"),
          nb::arg("guestIndices"),
          nb::arg("radius"));
    m.def("periodicCentroid",
          &site::periodicCentroid,
          "Centroid of a set of atoms with every atom unwrapped to its minimum image "
          "about the first.",
          nb::arg("yCloud"),
          nb::arg("atoms"));
    m.def("localTopologyKey",
          &topo::localKey,
          "Isomorphism-class key of the rooted neighbourhood of an atom within hops bonds "
          "(nauty certificate when linked, Weisfeiler-Lehman hash otherwise). colours is "
          "empty or one integer class per row; vertices of different colours never match.",
          nb::arg("rows"),
          nb::arg("atom"),
          nb::arg("hops") = 2,
          nb::arg("colours") = std::vector<int>{});
    m.def("topologyFingerprint",
          &topo::fingerprint,
          "Per-atom topology keys, their histogram, the primitive ring census and a "
          "label-independent frame key over neighbour rows by index; colours as in "
          "localTopologyKey.",
          nb::arg("rows"),
          nb::arg("hops") = 2,
          nb::arg("maxRingSize") = 7,
          nb::arg("colours") = std::vector<int>{});
    m.def(
        "parseSiteSpec",
        [](const std::string &spec) { return site::parseSiteSpec(spec); },
        "Parse '1=cationHead,2=anion[,family=ionicLiquid]'.",
        nb::arg("spec"));
    m.def("indicesOf",
          &site::indicesOf,
          "Cloud indices whose mapped kind matches (polar/apolar are unions).",
          nb::arg("yCloud"),
          nb::arg("table"),
          nb::arg("kind"));
    m.def("lammpsTypeOfKind",
          &site::lammpsTypeOfKind,
          "The unique LAMMPS type mapped to kind. Errors if not unique.",
          nb::arg("table"),
          nb::arg("kind"));
    m.attr("Kind") = m.attr("SiteKind");
    m.attr("Family") = m.attr("SiteFamily");
    nb::class_<site::DensityZ>(m, "DensityZ", "Cartesian number-density histogram.")
        .def_ro("z", &site::DensityZ::z)
        .def_ro("rho", &site::DensityZ::rho)
        .def_ro("type", &site::DensityZ::type);
    m.def(
        "densityZ",
        [](const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           int typeI,
           int nbin,
           int axis) { return site::densityZ(yCloud, typeI, nbin, axis); },
        "Type-resolved Cartesian number-density histogram.",
        nb::arg("yCloud"),
        nb::arg("typeI"),
        nb::arg("nbin"),
        nb::arg("axis"));
    m.def(
        "densityZ",
        [](const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
           const site::Table &table,
           site::Kind kind,
           int nbin,
           int axis) { return site::densityZ(yCloud, table, kind, nbin, axis); },
        "Site-kind-resolved Cartesian number-density histogram.",
        nb::arg("yCloud"),
        nb::arg("table"),
        nb::arg("kind"),
        nb::arg("nbin"),
        nb::arg("axis"));
    m.def("ionCloud",
          &site::ionCloud,
          "One COM vertex per ion molID, unwrapped with relDist. "
          "Output types are 1 (cationHead) and 2 (anion).",
          nb::arg("src"),
          nb::arg("table"));
    m.def(
        "ionCloud",
        [](const molSys::PointCloud<molSys::Point<double>, double> &src,
           int cationType,
           int anionType) {
            site::Table table;
            table.family = site::Family::ionicLiquid;
            table.typeToKind[cationType] = site::Kind::cationHead;
            table.typeToKind[anionType] = site::Kind::anion;
            return site::ionCloud(src, table);
        },
        "One COM vertex per ion molID from cation and anion LAMMPS types.",
        nb::arg("src"),
        nb::arg("cationType"),
        nb::arg("anionType"));
    m.def(
        "ionCloud",
        [](const molSys::PointCloud<molSys::Point<double>, double> &src,
           const std::unordered_map<int, site::Kind> &typeToKind) {
            site::Table table;
            table.typeToKind = typeToKind;
            return site::ionCloud(src, table);
        },
        "One COM vertex per ion molID from a type-to-Kind map.",
        nb::arg("src"),
        nb::arg("typeToKind"));
}
