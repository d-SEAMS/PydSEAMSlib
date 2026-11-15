
========
pydseams
========

.. raw:: html

   <p class="dseams-hero">
     <img class="dseams-hero-logo dseams-hero-logo--light"
          src="_static/logo/pydseams-logo-light.png"
          alt="pydseams"
          width="320"
          height="320"
          loading="eager" />
     <img class="dseams-hero-logo dseams-hero-logo--dark"
          src="_static/logo/pydseams-logo-dark.png"
          alt="pydseams"
          width="320"
          height="320"
          loading="eager" />
   </p>

.. grid:: 1 2 3 3
   :gutter: 2
   :padding: 1 1 0 0
   :class-container: sd-text-center

   .. grid-item-card:: Read
      :link: quickstart
      :link-type: doc
      :class-card: sd-shadow-sm

      Load a LAMMPS dump, XYZ, PDB, GRO, DCD, or eOn ``.con`` as a
      ``Frame``.

   .. grid-item-card:: Classify
      :link: tutorials/classify-ice
      :link-type: doc
      :class-card: sd-shadow-sm

      CHILL+ labels and HC / DDC cage scores on one frame.

   .. grid-item-card:: Bind
      :link: howto/ase
      :link-type: doc
      :class-card: sd-shadow-sm

      Round-trip ASE ``Atoms``. Optional solvis view of the same
      cloud.

About
-----

``pydseams`` is the Python package for the d-SEAMS engine. The compiled
module is ``yoda``. Helpers (``Frame``, ``read``, ASE, solvis) sit on that
surface. ``import pydseamslib`` still works.

Primary author: Ruhila S. The project started as PSF GSoC 2023
(``pyseams``; the PyPI name is ``pydseamslib``). Helpers stay in
Python so the C++ ABI does not grow a second I/O layer.

The engine and the ``seams`` CLI live in
`seams-core <https://github.com/d-SEAMS/seams-core>`_. Lua is
``dseams`` in `yodaStruct <https://github.com/d-SEAMS/yodaStruct>`_.
This site documents the Python ``Frame`` API. The header **Ecosystem**
menu jumps to the engine docs, the Lua front end, and the
neighbour backends (vesin cutoff lists, linkcell k-NN).

.. code:: python

    import pydseams as ds

    frame = ds.read("water.lammpstrj")
    print(frame.chill_plus())
    print(frame.cages())

Suite stack
-----------

Python helpers sit on the compiled ``yoda`` module. The same
``libyodaLib`` engine backs the ``seams`` CLI and the Lua ``dseams``
module.

.. mermaid::

   flowchart TB
     subgraph engine["seams-core"]
       LIB[libyodaLib]
     end
     subgraph compiled["pydseams.yoda"]
       YODA[yoda]
       ALIAS["_core / cyoda"]
     end
     subgraph helpers["Python helpers"]
       READ[ds.read]
       FRAME[Frame]
       ASE[from_ase / to_ase]
       SOL[to_solvis]
     end
     CLI["seams CLI"]
     LUA["yodaStruct / dseams"]
     LC["linkcell k-NN"]
     VE["vesin cutoff"]
     LC --> LIB
     VE --> LIB
     LIB --> YODA
     YODA --- ALIAS
     LIB --> CLI
     LIB --> LUA
     YODA --> READ
     YODA --> FRAME
     FRAME --> ASE
     FRAME --> SOL

Documentation structure
-----------------------

This documentation follows the `Diataxis <https://diataxis.fr/>`_
framework.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/index
   tutorials/classify-ice
   tutorials/dseams
   tutorials/bulk-ice
   tutorials/ice-and-clathrate
   tutorials/rings
   tutorials/order-parameters
   tutorials/confined-water
   tutorials/paper-figures

.. toctree::
   :maxdepth: 2
   :caption: How-to

   howto/index
   howto/install
   howto/ase
   howto/solvis
   howto/ase-integration
   howto/solvis-dseams
   howto/features
   howto/tools
   howto/faq
   howto/troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: Explanation

   explanation/index
   explanation/yoda-surface
   explanation/frame
   explanation/citation

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/index
   reference/python
   api

.. toctree::
   :maxdepth: 1
   :caption: Development

   history

Related projects
----------------

`seams-core <https://github.com/d-SEAMS/seams-core>`_
    ``libyodaLib`` and the ``seams`` CLI

`yodaStruct <https://github.com/d-SEAMS/yodaStruct>`_
    Lua / Fennel ``require("dseams")``

`linkcell <https://github.com/d-SEAMS/linkcell>`_
    periodic linked-cell k-nearest neighbours

vesin
    optional cutoff cell list (``neighListO``); brute force if not built

`d-SEAMS engine docs <https://docs.dseams.info>`_
    C++ API and theory

`dseams.info <https://dseams.info>`_
    project site

License
-------

MIT. Cite the 2020 d-SEAMS paper (DOI
`10.1021/acs.jcim.0c00031 <https://doi.org/10.1021/acs.jcim.0c00031>`_).
See `How to cite <explanation/citation.rst>`_.

#+begin\_export rst
Indices and tables
``================``

:ref:\`genindex\`
-----------------

:ref:\`modindex\`
-----------------

:ref:\`search\`
---------------

#+end\_export
