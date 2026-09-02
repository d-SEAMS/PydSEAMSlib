==============
Python surface
==============


Human lookup of ``Frame``, ``IceCounts``, ``CageScore``, ``DensityProfile``,
``ContactPairs``, ``DomainStats``, ``read``, the ASE helpers, and
``pydseams.features``. Napoleon autodoc of every signature is
`api.md <../../source/api.md>`_.

Frame
-----

``pydseams.frame.Frame`` is the one-frame handle. Load a LAMMPS dump,
an ASE ``Atoms``, or raw arrays, then call classification, density,
pairing, or domain methods. Classification does not write files.
Prefer ``pydseams.io.read``,
``Frame.from_ase``, or ``Frame.from_arrays`` over constructing ``Frame``
by filename.

``Trajectory`` is an alias of ``Frame``.

Constructors
~~~~~~~~~~~~

.. table::

    +--------------------------+-------------------------------------------------------+
    | name                     | role                                                  |
    +==========================+=======================================================+
    | ``Frame.from_file``      | LAMMPS dump; guesses type 2 then type 1               |
    +--------------------------+-------------------------------------------------------+
    | ``Frame.from_ase``       | ASE ``Atoms``; symbol, number, sequence, or ``None``  |
    +--------------------------+-------------------------------------------------------+
    | ``Frame.from_arrays``    | ``(N, 3)`` positions and three box lengths            |
    +--------------------------+-------------------------------------------------------+
    | ``Frame.from_xyz``       | XYZ via ``yoda.readXYZ``                              |
    +--------------------------+-------------------------------------------------------+
    | ``Frame.from_chemfiles`` | PDB / GRO / DCD when chemfiles is linked              |
    +--------------------------+-------------------------------------------------------+
    | ``Frame.from_con``       | eOn ``.con`` when readcon-core is linked              |
    +--------------------------+-------------------------------------------------------+

Package-level aliases: ``pydseams.from_ase``, ``from_arrays``,
``from_xyz``, ``from_chemfiles``, ``from_con``. Autodoc:
`api.md <../../source/api.md>`_.

Geometry and graphs
~~~~~~~~~~~~~~~~~~~

.. table::

    +-------------------+---------------------------------------------------------+
    | name              | role                                                    |
    +===================+=========================================================+
    | ``n_atoms``       | analysed particle count (``cloud.nop``)                 |
    +-------------------+---------------------------------------------------------+
    | ``box``           | three lengths or six LAMMPS restricted-triclinic values |
    +-------------------+---------------------------------------------------------+
    | ``positions``     | list of ``(x, y, z)``                                   |
    +-------------------+---------------------------------------------------------+
    | ``neighbor_list`` | cutoff neighbour list (``yoda.neighListO``)             |
    +-------------------+---------------------------------------------------------+
    | ``hbonds``        | hydrogen-bond neighbour list                            |
    +-------------------+---------------------------------------------------------+
    | ``rings``         | primitive rings up to size 6                            |
    +-------------------+---------------------------------------------------------+
    | ``load_frame``    | reload a later frame from the same file                 |
    +-------------------+---------------------------------------------------------+

``box`` is three lengths for an orthorhombic frame. Triclinic frames use
dump bound spans followed by ``xy``, ``xz``, and ``yz`` tilts.

``yoda.neighListPair``, ``yoda.SiteTable``, ``yoda.parseSiteSpec``,
``yoda.Kind`` / ``yoda.Family`` (aliases of ``SiteKind`` / ``SiteFamily``),
``yoda.ionCloud`` (``table``, ``cationType`` / ``anionType``, or
``typeToKind``), ``yoda.partialRdfHist``, ``yoda.coordinationNumber``,
``yoda.runningCN``, ``yoda.populateHbondsFromDonors``, and
``yoda.donatedHydrogenBond`` are bound. ``yoda.densityZ``,
``yoda.mutualNearestUnlike``, and ``yoda.largestDomain`` provide the raw
site-analysis primitives. ``Frame.rdf`` is ``(r, g)``;
``Frame.cn`` integrates that histogram; ``Frame.running_cn`` is the
running integral with ``rho_J = nJ / volume``; ``Frame.ion_cloud``
and ``Frame.hbonds_from_donors`` wrap the site and donor-H paths.
``Frame.density``, ``Frame.pairs``, and ``Frame.domain`` expose
type/site-resolved density, mutual cation--anion pairing, and the
largest connected site domain.

Classification
~~~~~~~~~~~~~~

.. table::

    +-------------------------+----------------------------------------------------+
    | name                    | role                                               |
    +=========================+====================================================+
    | ``chill_plus``          | CHILL+ ice labels; no files written                |
    +-------------------------+----------------------------------------------------+
    | ``chill``               | CHILL ice labels; no files written                 |
    +-------------------------+----------------------------------------------------+
    | ``classify_chill_plus`` | alias of ``chill_plus``                            |
    +-------------------------+----------------------------------------------------+
    | ``classify_chill``      | alias of ``chill``                                 |
    +-------------------------+----------------------------------------------------+
    | ``cages``               | HC / DDC membership; seeded hysteresis, ``k=4``    |
    +-------------------------+----------------------------------------------------+
    | ``seeded_affiliation``  | hysteresis flags; ``ring_adjacent`` completion     |
    +-------------------------+----------------------------------------------------+
    | ``ion_environment``     | first-shell ice class of every ion                 |
    +-------------------------+----------------------------------------------------+

``cages()`` with the default seeded path calls
``seeded_affiliation``. ``ring_adjacent`` fills the last vertex of a
six-ring whose other vertices carry a label. Pass ``seeded=0`` for
cutoff-graph affiliation on this frame's six-rings.

Export and descriptors
~~~~~~~~~~~~~~~~~~~~~~

.. table::

    +------------------------+--------------------------------------------------------+
    | name                   | role                                                   |
    +========================+========================================================+
    | ``to_ase``             | ASE ``Atoms``; needs ``pydseams[ase]``                 |
    +------------------------+--------------------------------------------------------+
    | ``to_solvis``          | solvis ``System``; needs ``pydseams[solvis]``          |
    +------------------------+--------------------------------------------------------+
    | ``steinhardt``         | local and neighbour-averaged ``ql``                    |
    +------------------------+--------------------------------------------------------+
    | ``soap``               | SOAP power spectrum                                    |
    +------------------------+--------------------------------------------------------+
    | ``voronoi_features``   | per-atom ``[q4, q6, q8]``                              |
    +------------------------+--------------------------------------------------------+
    | ``rdf``                | in-memory partial 3D RDF (``yoda.partialRdf``)         |
    +------------------------+--------------------------------------------------------+
    | ``cn``                 | site-site CN to a cutoff (``yoda.coordinationNumber``) |
    +------------------------+--------------------------------------------------------+
    | ``running_cn``         | running site-site CN (``yoda.runningCN``)              |
    +------------------------+--------------------------------------------------------+
    | ``density``            | Cartesian number density by particle type or site kind |
    +------------------------+--------------------------------------------------------+
    | ``ion_cloud``          | ion COM vertices (``yoda.ionCloud``)                   |
    +------------------------+--------------------------------------------------------+
    | ``pairs``              | mutual nearest unlike ion pairs                        |
    +------------------------+--------------------------------------------------------+
    | ``domain``             | largest cutoff-connected component of a site kind      |
    +------------------------+--------------------------------------------------------+
    | ``hbonds_from_donors`` | H-bond net from an explicit H-index list               |
    +------------------------+--------------------------------------------------------+

Prism, monolayer, and ``rdf_2d`` write engine output. ``Frame.rdf``
does not. See the Frame autodoc on `api.md <../../source/api.md>`_.

IceCounts
---------

``pydseams.frame.IceCounts`` is the histogram returned by
``chill_plus`` and ``chill``.

Keys are the ``AtomStateType`` names (``cubic``, ``hexagonal``, ``water``,
``interfacial``, ``clathrate``, ``interClathrate``, ``unclassified``,
``reCubic``, ``reHex``). Missing keys read as ``0`` via attribute access,
so ``counts.cubic`` and ``counts["hexagonal"]`` are equivalent.
``repr`` omits zero-count labels.

CageScore
---------

``pydseams.frame.CageScore`` is the per-atom score returned by
``cages``.

A molecule in an HC is ice Ih; a molecule in a DDC is ice Ic.

.. table::

    +-------------+-----------------------------------+
    | name        | role                              |
    +=============+===================================+
    | ``hc``      | per-atom hexagonal-cage flag      |
    +-------------+-----------------------------------+
    | ``ddc``     | per-atom double-diamond-cage flag |
    +-------------+-----------------------------------+
    | ``n_ih``    | number of atoms flagged HC        |
    +-------------+-----------------------------------+
    | ``n_ic``    | number of atoms flagged DDC       |
    +-------------+-----------------------------------+
    | ``n_water`` | number of atoms in neither cage   |
    +-------------+-----------------------------------+

Site-analysis records
---------------------

``DensityProfile``, ``ContactPairs``, and ``DomainStats`` are immutable
records returned by ``Frame.density``, ``Frame.pairs``, and ``Frame.domain``.

.. table::

    +--------------------+--------------------------------------------------------------+
    | record             | fields                                                       |
    +====================+==============================================================+
    | ``DensityProfile`` | ``centres``, ``rho``, ``axis``, ``atom_type``, ``site_kind`` |
    +--------------------+--------------------------------------------------------------+
    | ``ContactPairs``   | ``pairs``, ``count``, ``n_cation``, ``n_anion``              |
    +--------------------+--------------------------------------------------------------+
    | ``DomainStats``    | ``site_kind``, ``n``, ``largest``, ``percolation``           |
    +--------------------+--------------------------------------------------------------+

Use a site table when the analysis is chemistry-resolved:

.. code:: python

    frame = ds.read("ions.lammpstrj", all_atoms=True, atom_type=1)
    sites = ds.yoda.parseSiteSpec("1=polar,2=apolar")
    profile = frame.density(table=sites, kind=ds.yoda.Kind.polar)
    domain = frame.domain(sites, ds.yoda.Kind.polar)

    ions = ds.yoda.parseSiteSpec("1=cationHead,2=anion")
    pairs = frame.pairs(ions)

read
----

``pydseams.io.read`` is the suffix-dispatching loader.
``pydseams.available_readers`` reports which compiled readers this
build linked.

.. table::

    +-----------------------------------------------+--------------------------+
    | suffix                                        | constructor              |
    +===============================================+==========================+
    | ``.xyz``                                      | ``Frame.from_xyz``       |
    +-----------------------------------------------+--------------------------+
    | ``.con``                                      | ``Frame.from_con``       |
    +-----------------------------------------------+--------------------------+
    | ``.pdb``, ``.gro``, ``.dcd``                  | ``Frame.from_chemfiles`` |
    +-----------------------------------------------+--------------------------+
    | ``.lammpstrj``, ``.dump``, ``.lammps``, other | ``Frame.from_file``      |
    +-----------------------------------------------+--------------------------+

Common kwargs: ``frame`` (1-indexed), ``cutoff`` (Angstroms, default
3.5), ``bonded`` (``"auto"``, ``"hbond"``, ``"cutoff"``), ``atom_type``,
``region``, and ``all_atoms``. ``all_atoms=True`` retains every LAMMPS type
for mixed-site analyses; ``atom_type`` still selects the species used by
neighbour and ice methods. It cannot be combined with ``region``. For
LAMMPS readers, ``bonded="auto"`` uses hydrogen bonds when hydrogens are
available.

``region`` is ``(lo, hi)`` passed to ``yoda.readLammpsTrjreduced``.
``nop`` is the kept count. An axis with ``lo == hi`` is unconstrained,
so ``([0, 0, 0], [50, 0, 0])`` slices ``x`` only. ``yoda.readLammpsTrjO``
keeps every atom of the type and only sets ``inSlice``.

Builds without chemfiles or readcon raise ``RuntimeError`` when those
suffixes are used.

ASE
---

``pydseams.from_ase`` / ``Frame.from_ase`` build a ``Frame`` from an ASE
``Atoms``. ``Frame.to_ase`` is the inverse.

Install: ``pip install 'pydseamslib[ase]'``. The cell must be
nonsingular and periodic in all three directions. General cells are
rotated into LAMMPS restricted-triclinic form for analysis;
``to_ase`` restores the imported cell orientation and displacement.
``select`` is a symbol, an atomic number, a sequence of those, or
``None`` (every atom). Default ``select="O"``. A sequence such as
``("O", "Na", "Cl")`` keeps every listed species in the cloud and
analyses the first, so ions ride along for
``pydseams.features.ion_environment``.

For ASE input, ``bonded="auto"`` uses hydrogen bonds only when the
input contains ``H`` and the selected analysis cloud excludes ``H``.
Thus ``select=None`` uses cutoff bonding. Explicit ``bonded="hbond"``
requires a heavy-atom selection such as ``select="O"``.

``to_ase`` writes ``arrays["ice_type"]`` after CHILL and
``arrays["hc"]`` / ``arrays["ddc"]`` after ``cages()``. A frame loaded
from a LAMMPS dump (no ASE symbols) uses ``O`` as the fallback
species.

Implementation autodoc is the ASE helpers section of
`api.md <../../source/api.md>`_. Walkthrough:
`ASE how-to <../howto/ase.rst>`_.

``pydseams.to_solvis`` wraps the same ``Atoms``. Install:
``pip install 'pydseamslib[solvis]'``.

Features
--------

``pydseams.features`` builds per-frame kinetic vectors and
per-molecule states. ``IceFeaturizer`` walks a ``Trajectory`` (alias
of ``Frame``) and returns NumPy arrays. Autodoc:
`api.md <../../source/api.md>`_. Walkthrough:
`Features and ions <../howto/features.rst>`_.

.. table::

    +---------------------------+-----------------------------------------------------+
    | name                      | role                                                |
    +===========================+=====================================================+
    | ``IceFeaturizer``         | per-frame vector and per-molecule states            |
    +---------------------------+-----------------------------------------------------+
    | ``discretize_nmax``       | bin the largest-cluster series                      |
    +---------------------------+-----------------------------------------------------+
    | ``to_deeptime``           | time lagged independent component analysis (TICA)   |
    +---------------------------+-----------------------------------------------------+
    | ``to_pyemma_featurizer``  | register the vector with PyEMMA                     |
    +---------------------------+-----------------------------------------------------+
    | ``ion_environment``       | first-shell ice class of every ion                  |
    +---------------------------+-----------------------------------------------------+
    | ``ion_features``          | per-frame ion summary                               |
    +---------------------------+-----------------------------------------------------+

``IceFeaturizer`` turns on ring-adjacent completion
(``ring_adjacent``) of ``Frame.seeded_affiliation``. With
``ion_types`` the vector gains ``n_ion_ice``, ``n_ion_front``,
``n_ion_liquid`` and the mean shell ice fraction. ASE clouds that
must keep ions use ``Frame.from_ase(atoms, select=("O", "Na", "Cl"))``.

Package names
-------------

.. table::

    +-------------------------+-------------------------------------+
    | name                    | role                                |
    +=========================+=====================================+
    | ``pydseams``            | Python package                      |
    +-------------------------+-------------------------------------+
    | ``pydseamslib``         | compatibility alias of ``pydseams`` |
    +-------------------------+-------------------------------------+
    | ``pydseams.yoda``       | compiled registrations              |
    +-------------------------+-------------------------------------+
    | ``pydseams._core``      | alias of ``yoda``                   |
    +-------------------------+-------------------------------------+
    | ``pydseams.cyoda``      | alias of ``yoda``                   |
    +-------------------------+-------------------------------------+
    | ``pydseams.Trajectory`` | alias of ``Frame``                  |
    +-------------------------+-------------------------------------+
    | ``pydseams.features``   | kinetic vectors and ion environment |
    +-------------------------+-------------------------------------+

Requires Python 3.12+. Wheels are the CPython 3.12 limited ABI.

``yoda`` names (mocked in the Sphinx build) are tabulated on
`api.md <../../source/api.md>`_.
