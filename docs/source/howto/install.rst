================
Install pydseams
================


.. contents::


Problem
-------

You want ``import pydseams`` and ``import dseams`` in a Python environment,
optionally with ASE, solvis, or MDAnalysis integrations.

Requirements
------------

- Python 3.11 or later

- pip, pixi, or a source build toolchain (C++20 compiler, Meson, and
  Rust 1.88+ with ``cargo``)

Prebuilt wheels bundle the d-SEAMS C++ engine so no local compilation
is needed for standard installation.

pip
---

.. code:: bash

    pip install pydseamslib
    pip install 'pydseamslib[ase]'
    pip install 'pydseamslib[solvis]'
    pip install 'pydseamslib[mdanalysis]'

.. table::

    +-------------------------------------------+-------------------------------------------+
    | command                                   | provides                                  |
    +===========================================+===========================================+
    | ``pip install pydseamslib``               | ``pydseams``, ``pydseamslib``, ``dseams`` |
    +-------------------------------------------+-------------------------------------------+
    | ``pip install 'pydseamslib[ase]'``        | ``from_ase``, ``Frame.to_ase``            |
    +-------------------------------------------+-------------------------------------------+
    | ``pip install 'pydseamslib[solvis]'``     | ``to_solvis`` (includes ASE)              |
    +-------------------------------------------+-------------------------------------------+
    | ``pip install 'pydseamslib[mdanalysis]'`` | ``pydseams.adapters.IceStates``           |
    +-------------------------------------------+-------------------------------------------+

Without the corresponding extra, ``from_ase``, ``to_ase``, ``to_solvis``, and
``IceStates`` raise ``ImportError`` explaining which extra to install.

pixi
----

Add the PyPI package to a pixi project:

.. code:: bash

    pixi add --pypi pydseamslib

With extras:

.. code:: bash

    pixi add --pypi 'pydseamslib[ase]'
    pixi add --pypi 'pydseamslib[solvis]'

Verify
------

.. code:: bash

    python -c "import pydseams as ds, dseams; print(ds.__version__, dseams.__version__)"

``import pydseamslib`` re-exports ``pydseams``.

From source
-----------

The source build uses ``meson-python``. From a repository checkout:

.. code:: bash

    pip install -v '.[test,ase]'

With Nix:

.. code:: bash

    nix build
    nix develop

``nix build`` builds the package and ``nix develop`` opens a development
shell with ``pytest`` and ``hypothesis``.

Next steps
----------

`Quickstart <../quickstart.rst>`_
    classify your first frame

`Classify ice <../tutorials/classify-ice.rst>`_
    walk through a LAMMPS dump

`FAQ <faq.rst>`_
    names, extras, and formats
