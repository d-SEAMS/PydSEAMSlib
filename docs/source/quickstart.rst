==========
Quickstart
==========



Install
-------

#skip_lint_start

.. code:: bash

    pip install pydseamslib
    pip install 'pydseamslib[ase]'
    pip install 'pydseamslib[solvis]'

#skip_lint_end

``pydseams`` is the package. The ``[ase]`` extra pulls ASE for ``from_ase``
/ ``to_ase``. The ``[solvis]`` extra pulls ASE plus solvis-tools for
``to_solvis``.

With pixi, add the PyPI package:

.. code:: bash

    pixi add --pypi pydseamslib

Requires Python 3.11+. CPython 3.11 has a version-specific wheel. Python 3.12 and newer use
the CPython 3.12 limited ABI (one ``abi3`` wheel per platform). A wheel already links the engine; do not
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
Pass ``all_atoms`` on for mixed-type site analyses such as
``Frame.pairs`` and ``Frame.domain``.

``yoda`` is the compiled module. ``_core`` and ``cyoda`` are aliases of
``yoda``. Helpers (``Frame``, ``read``, ASE, solvis) stay in Python.

ASE ``Atoms``: ``ds.from_ase(atoms)``. See the `ASE
how-to <howto/ase.rst>`_. A worked mixed ice-water dump is the
`classify-ice tutorial <tutorials/classify-ice.rst>`_.

Seeded cages and mixed ASE selections
-------------------------------------

``Frame.seeded_affiliation`` is the hysteresis construction that
``cages`` with the seeded flag on calls. Ring-adjacent completion fills
the last vertex of a six-ring whose other vertices carry a label:

#skip_lint_start

.. code:: python

    print(frame.seeded_affiliation(ring_adjacent=True))
    print(frame.cages(seeded=True, ring_adjacent=True))

#skip_lint_end

``IceFeaturizer`` turns the same flag on by default.

Pass a sequence to ``from_ase`` when the ASE ``Atoms`` mix water and
salt. The listed species stay in the cloud:

#skip_lint_start

.. code:: python

    frame = ds.from_ase(atoms, select=("O", "Na", "Cl"))

#skip_lint_end

The first entry is the analysed water. The other species keep their
atomic numbers as ``c_type`` for ``ion_environment``. See the
`features how-to <howto/features.rst>`_.
