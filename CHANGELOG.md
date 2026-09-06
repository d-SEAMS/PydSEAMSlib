# Changelog

All notable changes to this project are documented in this file.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Unreleased notes live in [`changelog.d/`](changelog.d/) and are assembled
by [towncrier](https://towncrier.readthedocs.io/).

<!-- towncrier release notes start -->

## 2.10.0 (2026-09-06)

- Engine pinned at seams-core v2.10.0: F4 and Zeron q12, water-type CHILL, incomplete cages, TUM offload.


## 2.9.1 (2026-09-02)

- Engine pinned at seams-core v2.9.2: threaded cell-list neighbour lists.

## [2.9.0] - 2026-09-02

### Added

- `Frame.cages_by_signature` and `yoda.findBySignature`: closed polyhedra by ring-size census.
- Engine pinned at seams-core v2.9.0.
- `Frame.cages_by_signature` and `yoda.findBySignature`: closed polyhedra from a ring-size census (`sodalite`, `4:6,6:8`, named `hc` / `ddc` through the TUM finders). Engine wrap follows seams-core `884bed86`.
- `Frame.classify_topology` takes a sequence of libraries built at different hop counts and names each atom by the deepest that knows it; `LibraryMatch.depth` records which one did (`yoda.matchLibraries`).
- `Frame.guest_occupancy`, `yoda.guestOccupancy`, `yoda.periodicCentroid`: guests (methane, THF, ions) placed in enumerated cages.
- `KeyLibrary.coloured` and `FrameFingerprint.coloured` expose the colouring a library was built with.
- `pydseams.adapters.IceStates` (MDAnalysis analysis) and `pydseams.adapters.ice_states` (OVITO modifier function); a Colab notebook under `notebooks/`; PyPI keywords and classifiers.
- `IonEnvironment.members` lists each ion's shell; `Frame.hydration_shell_rings` and `yoda.shellRingCensus` count the rings through it by size. ASE hydrogen-bond adapters preserve an explicit ``mol-id`` array. Without molecule metadata, each hydrogen is assigned to its nearest selected atom under periodic minimum-image distance. Installed-package checks cover donor ownership and exported CHILL+/cage arrays. ASE all-atom imports use cutoff bonding. Hydrogen-bond topology requires a selection that excludes hydrogen; forcing ``bonded="hbond"`` on a mixed selection raises a clear error. `Frame(..., region=(lo, hi))` calls `readLammpsTrjreduced`, which drops atoms outside the box. `readLammpsTrjO` keeps every atom of the type and sets `inSlice`. An axis with `lo == hi` is unconstrained.

## [2.8.1] - 2026-09-02

### Changed

- Engine wrap pinned at seams-core v2.8.0; the 2.8.0 tag carried an empty revision.

## [2.8.0] - 2026-09-02

### Added

- Topology keys: `Frame.fingerprint(hops, max_ring_size, colour_types)`, `Frame.topology_library`, `Frame.classify_topology`; `yoda.topologyFingerprint`, `yoda.localTopologyKey`, `yoda.KeyLibrary`, `yoda.matchLibrary`.
- Ions: `Frame.ion_environment`, engine-backed `pydseams.features.ion_environment` and `yoda.ionEnvironment`; `Frame.cages(ring_adjacent=...)`.
- The CON reader is embedded (readcon-core v0.14.10 through the engine); docs for features, ions and sequence selection; autodoc names the last-vertex completion rule.
- Engine pinned at seams-core v2.8.0.

## [2.7.0] - 2026-09-02

### Added

- `Frame.seeded_affiliation(ring_adjacent=...)` exposes the engine's ring completion.
- `pydseams.features`: `IceFeaturizer` (per-frame feature vector and per-molecule states), `discretize_nmax`, `to_deeptime`, `to_pyemma_featurizer`.
- Ions: `ion_environment`, `ion_features`, `IceFeaturizer(ion_types=...)`; `Frame.from_ase(select=("O", "Na", "Cl"))` keeps listed species in the cloud as their atomic numbers.
- How-to `docs/source/howto/features.rst`. Engine pinned at seams-core v2.7.0; the test environment stays on CPython 3.12 or 3.13.

## [2.6.0] - 2026-08-17

### Changed

- `subprojects/seams-core.wrap` is `v2.6.0`.

## [2.5.1] - 2026-08-16

### Changed

- The PyPI distribution name is `pydseamslib`, the project that already has 0.0.1 / 0.0.2. `import pydseams` is unchanged. `pip install pydseamslib` is the install line.

## [2.5.0] - 2026-08-16

### Changed

- `subprojects/seams-core.wrap` is `v2.5.0`. The Python surface stays the 2.4.1 binds (`Frame.rdf`, `Frame.cn`, `Frame.running_cn`, `Frame.ion_cloud`, `Frame.hbonds_from_donors`). Ice-score `--family`, `seams cn --ions`, `seams pairs`, `seams domains`, and `seams density-z` live on the engine CLI.

## [2.4.1] - 2026-08-16

### Changed

- `subprojects/seams-core.wrap` stays `v2.4.0`. `yoda` binds `neighListPair`, `SiteTable` / `parseSiteSpec` / `ionCloud`, `partialRdfHist` / `coordinationNumber` / `runningCN` / `firstMinimumBin`, `populateHbondsFromDonors`, and `donatedHydrogenBond`. `yoda.Kind` / `yoda.Family` alias `SiteKind` / `SiteFamily`. `yoda.ionCloud` also accepts `(src, cationType, anionType)` or `(src, typeToKind)`. `Frame.cn`, `Frame.running_cn`, `Frame.ion_cloud`, and `Frame.hbonds_from_donors` call those. `Frame.rdf` stays `(r, g)`. `Frame.running_cn` uses `rho_J = nJ / volume` from the bound `PartialRdf`.

## [2.4.0] - 2026-08-16

### Changed

- `subprojects/seams-core.wrap` is `v2.4.0`. The engine adds `site::Table`, `rdf::coordinationNumber`, `populateHbondsFromDonors`, and `neighListPair`. pydseams still uses `neighList`, `populateHbonds`, and `rdf2Danalysis_AA`; those new symbols are not bound. Tilt dumps expose dump H on `cloud.box` (length >= 6). Orthorhombic frames stay three lengths. `pydseams.config` is installed with the wheel.

## [2.3.1] - 2026-08-16

### Changed

- `subprojects/seams-core.wrap` is `v2.3.1`. Remaining cutoff builders (`neighList`, `halfNeighList`, in-plane RDF) use vesin.

## [2.3.0] - 2026-08-16

### Changed

- `subprojects/seams-core.wrap` is `v2.3.0` (linkcell v0.3.0). Cutoff, frame, and *k* follow the engine twelve-factor table (`SEAMS_CONFIG` / `seams.env`, then the environment, then the argument). The docs mark is the hexagonal ice cage inside a `Frame`, as SVG.

## [2.2.5] - 2026-08-15

### Changed

- `subprojects/seams-core.wrap` is `v2.2.5` (linkcell v0.2.4). The Python signature is unchanged.

## [2.2.4] - 2026-08-15

### Changed

- `subprojects/seams-core.wrap` is `v2.2.4` (linked-cell k-nearest). A top-level `linkcell.wrap` is `v0.2.4` so the extension statically links the archive. `kNearestNeighbourList` keeps the same Python signature.

## [2.2.3]

### Changed

- Docs are orgmode plus Shibuya. Autodoc covers `Frame` and `yoda`.

## [2.2.2] - 2026-08-15

### Changed

- The compiled extension is `pydseams.yoda`. `_core` and `cyoda` remain aliases of that module.

## [2.2.1] - 2026-08-15

### Changed

- Flake-based Nix package for `pydseams`.

## [2.2.0] - 2026-08-15

### Changed

- The import is `pydseams`. `pydseamslib` remains a compatibility alias. This matches `dseams` (Lua) and the `seams` engine CLI.

## [2.1.0] - 2026-08-15

### Changed

- Python helpers sit on `_core` the way metatomic sits on its C surface. `ds.read` dispatches by suffix (LAMMPS, XYZ, chemfiles PDB/GRO/DCD, readcon `.con`). `from_chemfiles`, `from_con`, `from_xyz`, and `to_solvis` (`pydseams[solvis]`) are first-class.

## [2.0.1] - 2026-08-15

### Changed

- README documents `pip install` and `frame.to_ase()`.
- `to_ase()` writes `hc` / `ddc` from the last `cages()` result.
- `Frame()` defaults match `read()` (`atom_type` guessed, `bonded="auto"`).
- Sphinx release is 2.0.0. GSoC `pdm.lock` / `package.json` / pybind11 wrap are gone.

## [2.0.0] - 2026-08-15

### Changed

- nanobind rewrite against `seams-core` as a library wrap. Public surface is `Frame` / `read` / `from_ase`. Wheels are CPython 3.12 limited ABI (`abi3`). Primary author: Ruhila S. The Lua and Fennel CLI is [yodaStruct](https://github.com/d-SEAMS/yodaStruct).
