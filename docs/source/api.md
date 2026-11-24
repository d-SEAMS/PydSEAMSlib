# API reference

Napoleon autodoc of the Python helpers. Quick lookup of
{class}`~pydseams.frame.Frame`, {class}`~pydseams.frame.IceCounts`,
{class}`~pydseams.frame.CageScore`,
{class}`~pydseams.frame.DensityProfile`,
{class}`~pydseams.frame.ContactPairs`,
{class}`~pydseams.frame.DomainStats`, {func}`~pydseams.io.read`,
{mod}`pydseams.features`, and ASE is on the
[Python surface](reference/python) page.

The public constructors live on the {mod}`pydseams` package:
{func}`~pydseams.from_ase`, {func}`~pydseams.from_arrays`,
{func}`~pydseams.from_xyz`, {func}`~pydseams.from_chemfiles`,
{func}`~pydseams.from_con`, and {func}`~pydseams.to_solvis`.
{func}`~pydseams.has_chemfiles` and {func}`~pydseams.has_readcon`
report whether the PDB, GRO, and DCD reader and the eOn ``.con``
reader are linked.
{meth}`~pydseams.frame.Frame.to_ase` is the inverse of
{func}`~pydseams.from_ase`.

The compiled extension is the nanobind module `pydseams.yoda`.
`_core` and `cyoda` are compatibility aliases of `yoda`. `Trajectory`
is an alias of `Frame`.

Classification ({meth}`~pydseams.frame.Frame.chill_plus`,
{meth}`~pydseams.frame.Frame.chill`,
{meth}`~pydseams.frame.Frame.cages`) does not write files.
{meth}`~pydseams.frame.Frame.find_prisms`,
{meth}`~pydseams.frame.Frame.monolayer_rings`, and
{meth}`~pydseams.frame.Frame.rdf_2d` write engine output files.
{meth}`~pydseams.frame.Frame.rdf` computes in memory without writing
files.

```{eval-rst}
.. automodule:: pydseams
   :members:
   :imported-members:
   :exclude-members: Frame, IceCounts, CageScore, DensityProfile, ContactPairs, DomainStats, read, available_readers, yoda, _core, cyoda
```

## Frame

```{eval-rst}
.. automodule:: pydseams.frame
   :members:
   :exclude-members: read
   :show-inheritance:
   :member-order: bysource
```

## Features

Per-frame kinetic-model vectors and ion first-shell classes. See the
[features how-to](howto/features) for a walkthrough and the
[Python surface](reference/python) page for a summary table.

```{eval-rst}
.. automodule:: pydseams.features
   :members:
   :member-order: bysource
```

## I/O

```{eval-rst}
.. automodule:: pydseams.io
   :members:
   :member-order: bysource
```

## ASE helpers

Implementation of {func}`~pydseams.from_ase` and
{meth}`~pydseams.frame.Frame.to_ase`. Optional extra:
`pip install 'pydseamslib[ase]'`.

```{eval-rst}
.. automodule:: pydseams.aseio
   :members:
   :member-order: bysource
```

## solvis

Implementation of {func}`~pydseams.to_solvis` and
{meth}`~pydseams.frame.Frame.to_solvis`. Optional extra:
`pip install 'pydseamslib[solvis]'`.

```{eval-rst}
.. automodule:: pydseams.solvis
   :members:
   :member-order: bysource
```

## Compiled module (`pydseams.yoda`)

`pydseams.yoda` is the compiled nanobind extension module. The tables
below list every public symbol the bindings register. The published
wheels expose {{ yoda_names }} of them ({{ yoda_functions }} functions,
{{ yoda_types }} types, and the `Kind` / `Family` type aliases). The
counts use `src/bindings.cpp` and the build options for the wheels.
Four functions are compile-gated and absent from builds
without the corresponding optional engine dependency: `readChemfiles`
(chemfiles), `readCon` (readcon-core), `ira_match` and `sofi_point_group`
(libira).

Application code normally calls {class}`~pydseams.frame.Frame`,
{func}`~pydseams.io.read`, or {func}`~pydseams.from_ase`. Call
`pydseams.yoda` directly for low-level engine primitives.

### I/O

| name | role |
|------|------|
| `readXYZ` | Read atom coordinates from an XYZ file. |
| `readLammpsTrjreduced` | One LAMMPS frame, one type. With `isSlice`, drops atoms outside `(coordLow, coordHigh)`. `nop` is the kept count. An axis with `lo == hi` is unconstrained. |
| `readLammpsTrjO` | One LAMMPS frame, one type (the `O` is historical). With `isSlice`, sets `inSlice` and does not drop. `nop` is the type-filtered count. |
| `readLammpsTrj` | One LAMMPS frame, every atom type. |
| `readBonds` | Read bond connectivity from a formatted bond file. |
| `atomInSlice` | True when each component is in the closed interval, or that axis has `lo == hi`. |
| `readChemfiles` | Read any trajectory format supported by chemfiles (PDB, GRO, DCD, ...). Available when linked with chemfiles. |
| `readCon` | Read a `.con` file (eOn saddle-point search trajectories). Available when linked with readcon-core. |
| `writeDump` | Write a LAMMPS dump file for the current point cloud. |

### Neighbours and hydrogen bonds

| name | role |
|------|------|
| `clearNeighbourList` | Free memory for a neighbour list. |
| `getNewNeighbourListByIndex` | Build a neighbour list by index using a distance cutoff. |
| `kNearestNeighbourList` | Exact k-nearest bonded graph, union- or mutually-symmetrized. |
| `kNearestNeighbourPair` | Mutual and union k-nearest graphs from one candidate search. Returns `(mutual, union)`. |
| `periodicDistance` | Minimum Euclidean distance between cloud-indexed atoms in the periodic cell. |
| `shellSeparation` | Certificate pair (max k-th distance, min (k+1)-th distance) for the exact reduction of the k-nearest graph to a cutoff graph. |
| `halfNeighList` | Build a half neighbour list (each pair stored once) for one atom type. |
| `neighbourListByIndex` | Convert an atom-ID neighbour list to an index-based neighbour list. |
| `neighList` | Build a full neighbour list for two atom types within a cutoff. |
| `neighListO` | Build a full neighbour list for a single atom type within a cutoff. |
| `neighListPair` | I-J neighbour list (like-type reuses `neighListO`). |
| `SiteTable` / `parseSiteSpec` | Type-to-kind map. |
| `Kind` / `Family` | Aliases of `SiteKind` / `SiteFamily`. |
| `indicesOf` | Cloud indices whose mapped site kind matches, including polar/apolar unions. |
| `lammpsTypeOfKind` | Unique LAMMPS type mapped to a site kind. |
| `ionCloud` | One COM vertex per ion `molID`, unwrapped with `relDist`. Also `(src, cationType, anionType)` or `(src, typeToKind)`. |
| `mutualNearestUnlike` | Mutual nearest pairs between two particle types under periodic boundaries. |
| `partialRdfHist` | Partial 3D RDF as `PartialRdf` (r, g, count, volume, nI, nJ). |
| `runningCN` | Running site-site CN (`rhoJ` defaults to `nJ/volume`). |
| `coordinationNumber` | Site-site CN to `rMax` (`rhoJ` defaults to `nJ/volume`). |
| `populateHbondsFromDonors` | H-bond net from a flat list of `hCloud` indices. |
| `donatedHydrogenBond` | One donor-acceptor O-O-H test. |
| `createBondsFromCages` | Create bond connectivity from rings and cage information. |
| `getHbondDistanceOH` | Compute the O-H hydrogen bond distance between two atoms. |
| `populateHbonds` | Build the hydrogen-bond network from a trajectory and neighbour list. |
| `populateHbondsWithInputClouds` | Build hydrogen bonds from pre-loaded oxygen and hydrogen point clouds. |
| `trimBonds` | Remove duplicate bonds from a bond list. |

`Frame.rdf` wraps `partialRdf` and returns `(r, g)`. `Frame.cn`,
`Frame.running_cn`, `Frame.ion_cloud`, and `Frame.hbonds_from_donors`
wrap `coordinationNumber`, `runningCN`, `ionCloud`, and
`populateHbondsFromDonors`.

### Rings and cages

| name | role |
|------|------|
| `clearGraph` | Free memory for a graph object. |
| `countAllRingsFromIndex` | Find all possible rings (including non-shortest-path) up to `maxDepth`. |
| `ringNetwork` | Find all primitive (shortest-path) rings up to `maxDepth`. |
| `cageAffiliation` | Order-free per-ring cage classification: `(hc, ddc)` flag vectors. |
| `findBySignature` | Closed polyhedra matching a ring-size census or named table entry. |
| `seededCageAffiliation` | Seeded (hysteresis) per-atom cage flags: strict-graph seeds, permissive-graph completion. `ringAdjacentCompletion` fills the last vertex of a six-ring whose other vertices carry a label. |
| `RingUpdater` | Exact incremental primitive rings for a neighbour list. |
| `AffiliationUpdater` | Exact incremental per-ring cage classification for one frame. |
| `populateGraphFromIndices` | Create a graph object from an index-based neighbour list. |
| `populateGraphFromNListID` | Create a graph object from an atom-ID neighbour list and point cloud. |
| `removeNonSPrings` | Remove non-shortest-path rings using the Franzblau criterion. |
| `restoreEdgesFromIndices` | Restore graph edges from an index-based neighbour list. |

### Ring classification

| name | role |
|------|------|
| `assignPolygonType` | Assign atom types based on the ring size of n-membered rings. |
| `assignPrismType` | Assign atom types for atoms belonging to prism rings. |
| `clearRingList` | Free memory for a list of rings. |
| `compareRings` | Check whether two unordered rings contain the same elements. |
| `commonElementsInThreeRings` | Check whether three rings share at least one common element. |
| `deformedPrismTypes` | Get atom type values for deformed prisms. |
| `discardExtraTetragonBlocks` | Discard duplicate 4-membered ring pairs that are parallel in one dimension. |
| `findPrisms` | Identify which rings form prism blocks. |
| `findsCommonElements` | Return the common elements shared by two rings. |
| `findTripletInRing` | Search for a triplet of atoms within a ring. |
| `getSingleRingSize` | Extract rings of a specific size from a list of all rings. |
| `hasCommonElements` | Check whether two rings share any common elements. |
| `basalPrismConditions` | Test whether two rings satisfy strict basal prism conditions. |
| `relaxedPrismConditions` | Test whether two rings satisfy relaxed prism conditions (at least one bond). |
| `getEdgeMoleculesInRings` | Select edge molecules in rings that straddle the slice boundary. |
| `printSliceGetEdgeMoleculesInRings` | Select edge molecules in rings and write slice output files. |

### Topology writers

| name | role |
|------|------|
| `polygonRingAnalysis` | Classify rings in a quasi-2D monolayer and write output. |
| `bulkPolygonRingAnalysis` | Classify rings in a bulk system and write output. |
| `prismAnalysis` | Prism identification on rings up to `maxDepth`; writes output. |
| `rmAxialTranslations` | Remove axial translations from an ice nanotube for visualisation. |

### Topological unit matching

| name | role |
|------|------|
| `atomsFromCages` | Get atom indices belonging to cages in a given cluster. |
| `averageRMSDatom` | Average the per-atom RMSD over the number of shared cages. |
| `buildRefDDC` | Build a reference double-diamond cage from a template XYZ file. |
| `buildRefHC` | Build a reference hexagonal cage from a template XYZ file. |
| `clusterCages` | Cluster cages using Stillinger's algorithm and write XYZ output. |
| `shapeMatchDDC` | Shape-match a target double-diamond cage against a reference. |
| `shapeMatchHC` | Shape-match a target hexagonal cage against a reference. |
| `topoBulkCriteria` | Find HCs and DDCs in a bulk system using topological criteria. |
| `topoUnitMatchingBulk` | Run full topological unit matching for bulk water. |
| `updateRMSDatom` | Update per-atom RMSD from a cage shape-matching result. |

### Selection

| name | role |
|------|------|
| `getPointCloudOneAtomType` | Extract a point cloud containing only atoms of a given type. |
| `atomsInSingleSlice` | Mark atoms inside a rectangular volume slice. |
| `moleculesInSingleSlice` | Mark whole molecules as in-slice if any atom is inside the region. |
| `setAtomsWithSameMolID` | Set the `inSlice` flag for all atoms sharing a given molecule ID. |

### CHILL / CHILL+

| name | role |
|------|------|
| `BondClassifier` | Bond-classification rule set (staggered / eclipsed windows). |
| `chillRule` | The CHILL water rule set. |
| `chillPlusRule` | The CHILL+ water rule set. |
| `bondClassifier` | Look up a registered bond-classification rule set by name. |
| `registerBondClassifier` | Register (or replace) a named bond-classification rule set. |
| `bondClassifierNames` | Names of every registered bond-classification rule set. |
| `classifyBonds` | Compute and classify bond correlations under an arbitrary rule set. |
| `getCorrelPlus` | Compute CHILL+ bond-order correlations and classify bond types. |
| `getIceTypePlusNoPrint` | Classify each atom's ice type using CHILL+. Does not write a file. |
| `getIceTypePlus` | Classify each atom's ice type using CHILL+ and write to file. |
| `getCorrel` | Compute CHILL bond-order correlations and classify bond types. |
| `getIceTypeNoPrint` | Classify each atom's ice type using CHILL. Does not write a file. |
| `getIceType` | Classify each atom's ice type using CHILL and write to file. |
| `getq6` | Compute the q6 bond order parameter for all atoms. |
| `reclassifyWater` | Reclassify water molecules using averaged q6 and q3 parameters. |
| `printIceType` | Print the ice type classification for the current frame. |

### Descriptors

| name | role |
|------|------|
| `steinhardtQl` | Local and neighbour-averaged Steinhardt parameters of degree `orderL`. |
| `steinhardtQlVoronoi` | Voronoi facet-area weighted Steinhardt parameters. |
| `classifyTemplates` | IRA/Horn overlay onto FCC, HCP, BCC, and SC neighbour shells. |
| `soapSpectrum` | SOAP power spectrum of one particle. |
| `soapSpectrumAll` | SOAP power spectrum of every particle. |
| `voronoiFeature` | Per-atom `[q4, q6, q8]` from the Voronoi-weighted Steinhardt path. |
| `voronoiFeatures` | `[q4, q6, q8]` for every particle from one Voronoi pass per order. |
| `voronoiFacetWeights` | Voronoi facet neighbours and area weights for every particle. |
| `LinearClassifier` | Linear classifier used by `Frame.fit_classifier`. |
| `ira_available` | True when this build linked libira (IRA/SOFI). |
| `ira_match` | IRA overlay of two `n x 3` point sets. Available when linked with libira. |
| `sofi_point_group` | SOFI point group of an `n x 3` cloud. Available when linked with libira. |
| `lookupTableQ4Vec` | Lookup table for Q4 (`m=0` to `m=8`). |
| `lookupTableQ4` | Lookup table for Q4 at a single `m` (`m=0` to `m=8`). |
| `lookupTableQ8Vec` | Lookup table for Q8 (`m=0` to `m=16`). |
| `lookupTableQ8` | Lookup table for Q8 at a single `m` (`m=0` to `m=16`). |

### Clustering and RDF

| name | role |
|------|------|
| `clusterAnalysis` | Cluster ice-like particles and return the largest ice cluster. |
| `recenterClusterCloud` | Recenter a cluster point cloud for visualisation. |
| `rdf2Danalysis_AA` | 2D radial distribution function for identical atom types. |
| `partialRdf` | Partial 3D RDF `g_IJ(r)`. Returns `(r, g)`. |
| `firstMinimumBin` | Bin index of the first minimum of g after the first peak, or -1. |
| `runningCN` | Running CN `4 pi rho_J int s^2 g ds`. Returns a list. |
| `densityZ` | Cartesian number-density histogram by particle type or mapped site kind and axis. |
| `largestDomain` | Size and fraction of the largest connected component in a selected subset. |

### Sites, guests, and local topology

| name | role |
|------|------|
| `ionEnvironment` | Class each ion by its first water shell against an ice flag. |
| `shellRingCensus` | Census, by size, of rings with at least one vertex in an ion's first shell. |
| `localTopologyKey` | Isomorphism-class key of one atom's rooted neighbourhood within `hops` bonds. |
| `topologyFingerprint` | Per-atom topology keys, their histogram, the primitive ring census, and a frame key. |
| `addToLibrary` | Add every distinct key of a fingerprint under a label. |
| `writeLibrary` | Text form of a key library. |
| `readLibrary` | Key library from its text form. |
| `matchLibrary` | Name every atom of a fingerprint by one key library. |
| `matchLibraries` | Name every atom by key libraries at several hop counts. |
| `guestOccupancy` | Assign each guest to the nearest cage centroid within a radius. |
| `guestOccupancyInside` | Assign each guest to a cage whose faces contain it. |
| `periodicCentroid` | Centroid of a set of atoms, each unwrapped to its minimum image about the first. |

### Types

| name | role |
|------|------|
| `PartialRdf` | Histogram from `partialRdfHist`: `r`, `g`, `count`, `volume`, `nI`, `nJ`. |
| `PointCloudDouble` | Collection of points for a single frame, with box dimensions. |
| `PointDouble` | Per-particle data: coordinates, type, molecule ID, ice classification. |
| `SteinhardtQl` | Per-particle Steinhardt `ql` and neighbour-averaged `qlBar`. |
| `AtomStateType` | Per-atom ice phase classification from CHILL / CHILL+ / q6. |
| `BondType` | Bond classification: staggered, eclipsed, or `out_of_range`. |
| `Result` | Bond correlation result: `classifier` (bond type) and `c_value`. |
| `CrystalKind` | Crystal template kind: `other`, `sc`, `fcc`, `hcp`, `bcc`. |
| `TemplateHit` | IRA/Horn overlay hit: `kind`, `rmsd`, `name`. |
| `VoronoiWeights` | Facet-sharing neighbours and facet-area weights of one particle. |
| `DensityZ` | Number-density result: bin centres `z`, values `rho`, and selected `type`. |
| `Domain` | Connected-domain result: subset size, largest component, and percolation fraction. |
| `SiteKind` | Per-atom site chemistry. Type 1 is not a kind. |
| `SiteFamily` | Input family. It is not inferred from types. |
| `IonState` | An ion by its first water shell: `liquid`, `front`, or `ice`. |
| `IonEnvironment` | Ions read against a per-atom ice assignment. |
| `LocalKey` | Key of one atom's rooted bonded neighbourhood. |
| `FrameFingerprint` | Local keys of every atom, their histogram, the ring census, and one frame key. |
| `KeyLibrary` | Local keys of reference structures with their labels. |
| `LibraryMatch` | Atoms named by a key library. |
| `GuestOccupancy` | Guests placed in enumerated cages by their periodic centroids. |
