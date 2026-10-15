================
Install pydseams
================


.. contents::


Problem
-------

You want ``import pydseams`` in a working interpreter, with optional
ASE or solvis extras.

Requirements
------------

- Python 3.11+

- pip, pixi, or a meson-python checkout

CPython 3.11 has a version-specific wheel. Python 3.12 and newer use
the CPython 3.12 limited ABI (one ``abi3`` wheel per
platform). A wheel already links the d-SEAMS engine. Do not compile
``yoda`` to install the package.

pip
---

.. code:: bash

    pip install pydseamslib
    pip install 'pydseamslib[ase]'
    pip install 'pydseamslib[solvis]'

.. table::

    +---------------------------------------+--------------------------------+
    | extra                                 | provides                       |
    +=======================================+================================+
    | ``pip install pydseamslib``           | ``Frame``, ``read``, ``yoda``  |
    +---------------------------------------+--------------------------------+
    | ``pip install 'pydseamslib[ase]'``    | ``from_ase``, ``Frame.to_ase`` |
    +---------------------------------------+--------------------------------+
    | ``pip install 'pydseamslib[solvis]'`` | ``to_solvis`` (pulls ASE)      |
    +---------------------------------------+--------------------------------+

The ``[solvis]`` extra installs ASE as well. Without an extra,
``from_ase`` / ``to_ase`` / ``to_solvis`` raise ``ImportError`` and name
the matching ``pip install`` command.

pixi
----

Add the PyPI package to a pixi project:

.. code:: bash

    pixi add --pypi pydseamslib

Extras:

.. code:: bash

    pixi add --pypi 'pydseamslib[ase]'
    pixi add --pypi 'pydseamslib[solvis]'

Verify
------

.. code:: bash

    python -c "import pydseams as ds; print(ds.__version__); print(ds.yoda is ds._core)"

The first line prints the installed version. The second line prints
``True``: ``_core`` and ``cyoda`` are aliases of ``yoda``.

``import pydseamslib`` is a compatibility alias of ``pydseams``.

From a checkout
---------------

The package build is meson-python. From the repository root:

.. code:: bash

    pip install -v '.[test,ase]'
    pip install -e .

Nix, at the repository root:

.. code:: bash

    nix build
    nix develop

``nix build`` produces the ``pydseams`` package. ``nix develop`` is the
dev shell (pytest, hypothesis). There is no CMake install.

The published wheel is the install path for the docs site and for
application code.

Next steps
----------

`Quickstart <../quickstart.rst>`_
    classify one frame

`Classify ice <../tutorials/classify-ice.rst>`_
    fixture dump

`FAQ <faq.rst>`_
    names and extras
