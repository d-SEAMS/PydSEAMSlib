================
The yoda surface
================


The import package is ``pydseams``. The compiled extension module
inside it is ``yoda`` (``pydseams.yoda``). That separation is intentional.

The ``yoda`` name
-----------------

d-SEAMS originally shipped in
`yodaStruct <https://github.com/d-SEAMS/yodaStruct>`_ (Goswami,
Goswami, Singh, *Journal of Chemical Information and Modeling* 2020,
DOI `10.1021/acs.jcim.0c00031 <https://doi.org/10.1021/acs.jcim.0c00031>`_).
The underlying C++ library is ``libyodaLib`` in
`seams-core <https://github.com/d-SEAMS/seams-core>`_, and the Lua
and Fennel front end lives in yodaStruct (``require("dseams")``).

The nanobind module is named ``yoda`` because it wraps that compiled C++
surface directly. ``pydseams`` is the Python package around it;
``import pydseamslib`` is a compatibility alias that re-exports
``pydseams``.

Why helpers stay in Python
--------------------------

``yoda`` exposes the compiled C++ primitives: ``PointCloudDouble``, file
readers, neighbour lists (``neighList``, ``neighListO``, ``neighListPair``),
CHILL and CHILL+, rings, cages, Steinhardt order parameters, SOAP,
``SiteTable``, ``partialRdfHist``, ``coordinationNumber``, ``runningCN``, and
``populateHbondsFromDonors``.

High-level ergonomics (``Frame``, ``read``, ``from_ase``, ``from_arrays``,
``from_xyz``, ``from_chemfiles``, ``from_con``, and ``to_solvis``) remain in
Python on top of ``yoda``. Keeping suffix dispatch and third-party
adapters in Python lets optional integrations (``pydseamslib[ase]``,
``pydseamslib[solvis]``, ``pydseamslib[mdanalysis]``) stay decoupled from
the compiled C++ ABI.

``ds.read`` dispatches on the file suffix and calls the matching
``Frame`` constructor. ``Frame(..., region=(lo, hi))`` delegates to
``yoda.readLammpsTrjreduced``, which drops atoms outside the box;
``yoda.readLammpsTrjO`` keeps every atom of the type and marks ``inSlice``.
Both remain directly callable when you need the raw point cloud.

Module aliases
--------------

.. code:: python

    import pydseams as ds

    assert ds._core is ds.yoda
    assert ds.cyoda is ds.yoda

.. table::

    +-------------------------+---------------------------------+
    | name                    | role                            |
    +=========================+=================================+
    | ``pydseams.yoda``       | compiled module (canonical)     |
    +-------------------------+---------------------------------+
    | ``pydseams._core``      | compatibility alias of ``yoda`` |
    +-------------------------+---------------------------------+
    | ``pydseams.cyoda``      | compatibility alias of ``yoda`` |
    +-------------------------+---------------------------------+
    | ``pydseamslib``         | package alias of ``pydseams``   |
    +-------------------------+---------------------------------+
    | ``pydseams.Trajectory`` | alias of ``Frame``              |
    +-------------------------+---------------------------------+

New code should import ``pydseams`` and access ``pydseams.yoda`` when
calling the compiled surface directly.
