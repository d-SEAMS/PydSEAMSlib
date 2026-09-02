==========
Quickstart
==========



Install
-------

.. code:: bash

    pip install pydseamslib
    pip install 'pydseamslib[ase]'
    pip install 'pydseamslib[solvis]'

``pydseams`` is the package. The ``[ase]`` extra pulls ASE for ``from_ase``
/ ``to_ase``. The ``[solvis]`` extra pulls ASE plus solvis-tools for
``to_solvis``.

With pixi, add the PyPI package:

.. code:: bash

    pixi add --pypi pydseamslib

Requires Python 3.12+. Wheels are the CPython 3.12 limited ABI (one
``abi3`` wheel per platform). A wheel already links the engine; do not
compile ``yoda`` to use the package.

``import pydseamslib`` is a compatibility alias of ``pydseams``.

From a PydSEAMSlib checkout:

.. code:: bash

    nix build
    nix develop

``nix build`` produces the ``pydseams`` package. ``nix develop`` is the
dev shell (pytest, hypothesis).

The full extra table is in the `install
how-to <howto/install.rst>`_.

Classify a frame
----------------

.. code:: python

    import pydseams as ds

    frame = ds.read("water.lammpstrj")   # also .xyz, .pdb, .gro, .dcd, .con
    print(frame.chill_plus())
    print(frame.cages())

``ds.read`` picks the engine reader from the suffix and returns a
``Frame``. ``chill_plus`` and ``cages`` do not write files.
Pass ``all_atoms=True`` for mixed-type site analyses such as
``Frame.pairs`` and ``Frame.domain``.

``yoda`` is the compiled module. ``_core`` and ``cyoda`` are aliases of
``yoda``. Helpers (``Frame``, ``read``, ASE, solvis) stay in Python.

Those ASE helpers take ``ds.from_ase(atoms)``. A sequence such as
``select=("O", "Na", "Cl")`` keeps listed species in the cloud and
analyses the first, so ions ride along. ``cages()`` calls
``seeded_affiliation``; ``ring_adjacent`` fills the last vertex of
a six-ring whose other vertices carry a label. See the `ASE
how-to <howto/ase.rst>`_. A worked mixed ice-water dump is the
`classify-ice tutorial <tutorials/classify-ice.rst>`_. The kinetic
vector path is the `features how-to <howto/features.rst>`_.
