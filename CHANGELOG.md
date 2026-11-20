# Changelog

All notable changes to this project are documented in this file.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Unreleased notes live in [`changelog.d/`](changelog.d/) and are assembled
by [towncrier](https://towncrier.readthedocs.io/).

<!-- towncrier release notes start -->

## [2.14.1] - 2026-11-19

### Added

- README and the docs landing page link the d-SEAMS book (Python, Python Frame) and the walkthrough video.

### Changed

- The documentation site moves to https://python.dseams.info/; the Lua links point at https://lua.dseams.info/.
- With seams-core 2.14.1, `Frame.rdf_2d` returns zero-valued bins for a frame whose atoms are all coplanar, matching the engine's empty-shell convention.

### Fixed

- A `dseams.userframe` frame with no atoms returns a `(0, 3)` coordinate array, a zero in-plane g(r), and zero prism rows. The coordinate property used to raise a null-pointer error, and the in-plane RDF and prism census crashed inside the engine.
- An all-liquid frame has no brine pockets or reservoir waters.
- The `libseams_lang` mutual and union graph paths discard k-nearest candidates beyond the requested cutoff plus 1.5 instead of using that distance only as a search-cell hint.
- Malformed MOL2 and mmCIF files that chemfiles refuses now report a read error through `dseams` instead of returning a frame with zero atoms.
- Links between docs pages written as Org file links resolve to the built HTML pages instead of missing .rst and .md files.
- The dseams reference matches the live signatures. `dseams.cages` takes `graph="seeded"`. `dseams.signature_cages` takes `graph="cutoff"` and defaults `spec` to `"512"`. `dseams.symmetry` takes `shell_rmsd` (default `0.12`).


## [2.14.0] - 2026-11-17

### Added

- A contributing guide at the repository root, `llms.txt` (also served at the docs site root), and a README section on developing against seams-core and the release order.
- The Python tutorial and how-to pages, examples and ctypes checks move here from seams-core: bulk ice, ice and clathrate, rings, order parameters, confined water, paper figures, ASE integration, solvis, the `dseams` autodoc and language-task pages, `lang_tasks_python.py` and `science_tutorials.py`. `tests/dseams/test_doc_snippets.py` runs their Python blocks against the installed wheel.
- `Frame.guest_occupancy_both(signature, guest_types, radius)` returns the radius and ray-parity occupancy of the same guests in the same signature cages, and `yoda.guestOccupancyInside` binds the engine's inside assignment.
- `import dseams` ships in the `pydseamslib` wheel: the nanobind module over the seams-core 2.14 C API (`chill_plus`, `q4q6`, `neighbours`, `rings`, `dseams.Frame`), `cages`, `signature_cages` and `symmetry` through `libseams_lang`, and `dseams.userframe`. The sources and tests move here from seams-core with their history; `libseams` and `libseams_lang` are inside the wheel.
- `src/seams` is the ctypes package over `libseams`, with the Python examples under `examples/python` and `examples/ctypes`. It is not in the wheel; `meson test` runs it against the subproject build.

### Changed

- The engine is seams-core 2.14.0. `FrameFingerprint.colored` and `KeyLibrary.colored` join the `coloured` names.
- The engine wrap follows seams-core tag `v2.14.0`. That release builds the C API examples and fuzzers only with `with_tests`, adds `seams_rings_chordless`, and requires linkcell `v0.4.2`.
- linkcell is v0.4.2, the release seams-core 2.14.0 pins. Meson takes a dependency both projects wrap from this project's wrap, so the v0.4.0 pin stopped the engine at its `linkcell >= 0.4.1` requirement once the engine moved past 0.4.0.
- The package metadata states the licence as the SPDX expression `MIT` with `LICENSE` as its licence file (core metadata 2.4, PEP 639), replacing the free-text licence and the licence classifier.
- This repository builds and publishes `pydseamslib` 2.14.0; seams-core no longer does. The sdist carries every subproject, the wheels are built from it on Linux (x86_64, aarch64) and macOS (arm64, x86_64), and PyPI trusted publishing runs from `build_wheels.yml` on a `v*` tag of d-SEAMS/PydSEAMSlib.
- `pydseams` CHILL+ counts on cutoff neighbour lists change with the engine: seams-core 2.14.0 fixed the shell count, so `Frame.chill_plus()` on `tests/data/exampleTraj.lammpstrj` (oxygen, type 2) reports 9 interClathrate and 241 water where 2.11 reported 12 and 238.

### Fixed

- On macOS the `core_counts`, `lang_tasks_python`, `minimal_py`, `python_api` and `science_tutorials` examples, and the tests that run them, load `libseams.dylib` and `libseams_lang.dylib` instead of the Linux `.so` names, so the documentation checks pass against the macOS wheels.
- The Nix flake builds seams-core tag `v2.14.0`, the revision `subprojects/seams-core.wrap` pins, with linkcell 0.4.2 and readcon-core 0.16.2 as the engine requires.
- The citation page lists the d-SEAMS paper in a bibliography so its `:cite:` resolves, and the ASE cross-references load from `https://docs.ase-lib.org/`, where the ASE inventory moved.
- The documentation is regenerated from its Org sources: the stray `#skip_lint_start` and `#skip_lint_end` lines are gone from the quickstart, ASE and features pages, the Python surface tables render instead of failing on an over-wide cell, and the hand edits (the CPython 3.11 wheel, `hydration_shell_rings`, `guest_occupancy`, `classify_topology` and the tools page) now live in the Org sources.
- The published API reference documents `dseams` and `pydseams.features` again; the documentation build installs NumPy, which autodoc needs to import them.
- The README now pins the released package in its install examples and describes the import names and readers that the 2.14.0 wheel actually contains.
- The macOS Nix build gives minimage, linkcell and readcon-core a store install name, so importing `pydseams.yoda` loads those libraries from the Nix store.
- The macOS Nix build passes Cargo's rustc probe. The Darwin install-name flag is an argument to rustc, after cargo-auditable.
- In exclude mode, `run_unittest.py` parses each selected module to find pytest-style tests, so a `test_*` function runs even when its module never imports pytest. A meson self-test covers the routing cases.
- `meson test python_dseams_no_checkout` runs in both engine modes and covers every `tests/dseams/test_*.py` module except `test_seams_package.py` and `test_doc_snippets.py`. The checkout-dependent pair stays on the subproject-only `python_dseams` test.
- `meson test` no longer kills `python_dseams` and `python_dseams_no_checkout` with signal 11 or 7 at random. Meson starts the two at once. Each runs `tests/dseams/run_unittest.py`, and that script copied `dseams._core` and the Python sources over the same files in `python-stage`. Rewriting a shared object in place corrupts any process that has it mapped, and `MALLOC_PERTURB_` played no part. The script now writes each file under a temporary name and renames it into place. The new `python_dseams_stage_concurrency` test starts twelve runs at once for ten rounds and checks that a published file gets a new inode.


## [2.11.1] - 2026-10-15

### Fixed

- Build CPython 3.11 wheels with its native API while retaining stable-ABI wheels for Python 3.12 and newer.


## [2.11.0] - 2026-10-09

### Added

- `pydseams.band.band_rows` names the cutoff pairs that enter or leave between frames, with the shortest same-type distance and the primitive-ring census.

### Changed

- CHANGELOG.md follows Keep a Changelog. Unreleased notes are
  changelog.d fragments assembled by towncrier.
- Engine wrap and flake lock follow seams-core 67628524 (v2.10.0 plus ice XXI).
- Engine wrap follows seams-core 91235ec1 (CHILL+ device-safe target loops).
- Engine wrap follows seams-core bdc9b74 (RDF MIC-once and ice XXI fingerprint).
- Engine wrap follows seams-core eb1a4985 (j815 leftover batch).
- Require minimage and linkcell in the bindings and use seams-core 2.11.0 for periodic distances, threaded neighbour lists, and topology keys.
- Wrap seams-core fa1c8a7397d0f9b020da07489e4211fa501f8d6d (IRA residual maps target onto ref).

### Fixed

- Pin the top-level geometry wraps to the engine versions so Meson resolves the same libraries.

- Band distances respect each frame's periodicity and use Euclidean nearest images for skew cells. Band comparisons require consistent atom IDs, types, and order.
- Build the Rust geometry libraries inside the wheel containers.
- `Frame.from_xyz` puts every element into the cutoff graph. Each symbol gets its own `c_type`. A nitrogen pair inside the cutoff stays bonded. `periodic=False` keeps that graph in free space, so a molecule is not wrapped across its bounding span.
- `pydseams.__version__` is 2.10.0, matching the package. `kNearestNeighbourList` binds both the type-I and type-set overloads so the wrap builds against seams-core after ice XXI.


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

- nanobind rewrite against `seams-core` as a library wrap. Public surface is `Frame` / `read` / `from_ase`. Wheels are CPython 3.12 limited ABI (`abi3`). Primary author: Ruhila S (https://orcid.org/0000-0002-5443-9356). The Lua and Fennel CLI is [yodaStruct](https://github.com/d-SEAMS/yodaStruct).
