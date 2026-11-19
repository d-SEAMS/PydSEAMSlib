==========================
Frequently asked questions
==========================


.. contents::


Names
-----

Why is the package ``pydseams`` and the compiled module ``yoda``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``pydseams`` is the Python package, and ``yoda`` (``pydseams.yoda``) is the
compiled nanobind extension inside it, named after the underlying
``libyodaLib`` engine in seams-core. High-level helpers (``Frame``, ``read``,
ASE, solvis) sit in Python on top of ``yoda``.

See `The yoda surface <../explanation/yoda-surface.rst>`_.

What is ``pydseamslib``?
~~~~~~~~~~~~~~~~~~~~~~~~

The PyPI distribution name (``pip install pydseamslib``).
``import pydseamslib`` is an alias of ``pydseams``. New code imports
``pydseams`` (or ``dseams`` for the array API).

What are ``_core``, ``cyoda``, and ``Trajectory``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``_core`` and ``cyoda`` are compatibility aliases of ``pydseams.yoda``
(``ds._core is ds.yoda`` and ``ds.cyoda is ds.yoda``). ``pydseams.Trajectory``
is an alias of ``pydseams.Frame``.

Do I need to compile ``yoda`` to use the package?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No. ``pip install pydseamslib`` installs a wheel that bundles the
compiled d-SEAMS engine. Compiling from source is only needed when
developing the bindings or building with optional engine libraries
such as chemfiles.

Usage
-----

Which formats does ``ds.read`` accept?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``ds.read`` dispatches by file suffix:

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

``available_readers()`` reports which optional C++ readers your build
includes (``xyz``, ``chemfiles``, ``readcon``). ``has_chemfiles()`` and
``has_readcon()`` are the same two flags. ``lammps`` is always enabled.

Does classification write files?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``chill_plus``, ``chill``, and ``cages`` do not write files. ``find_prisms``,
``monolayer_rings``, and ``rdf_2d`` write engine output files to the
directory passed to them.

Can ASE cells be triclinic?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes. ``from_ase`` rotates any nonsingular, fully periodic ASE cell into
the engine's LAMMPS restricted-triclinic representation, and ``to_ase``
restores the original cell orientation, positions, and displacement.
Partially periodic and singular cells are rejected.

How do I keep every atom, not just oxygen?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    frame = ds.from_ase(atoms, select=None)

The default is ``select="O"``. A symbol or atomic number selects that
species. A sequence such as ``("O", "Na", "Cl")`` keeps all listed
species and analyses the first. Mixed selections containing ``H`` use
cutoff bonding; hydrogen-bond topology requires a heavy-atom selection
such as ``select="O"``. For LAMMPS files, pass ``all_atoms=True`` to
``ds.read``.

How do I view a classified frame?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use solvis: ``pip install 'pydseamslib[solvis]'``, then
``frame.to_solvis()``. See `View a frame in solvis <solvis.rst>`_.

How do I cite d-SEAMS?
~~~~~~~~~~~~~~~~~~~~~~

Cite the 2020 d-SEAMS paper. See
`How to cite <../explanation/citation.rst>`_.

Compatibility
-------------

Which Python versions are supported?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python 3.11 and later (Linux and macOS). Free-threaded CPython is not
supported.

How is ``pydseams`` related to ``seams`` and Lua ``dseams``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

They share the same C++ engine (``libyodaLib``). ``seams`` is the CLI in
seams-core, ``require("dseams")`` is the Lua / Fennel interface in
yodaStruct, and ``pydseams`` / ``dseams`` are the Python bindings.

See also
--------

- `Install pydseams <install.rst>`_

- `Troubleshooting <troubleshooting.rst>`_

- `Classify ice <../tutorials/classify-ice.rst>`_
