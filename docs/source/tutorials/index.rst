=========
Tutorials
=========


Tutorials are learning-oriented. Follow them in order to learn the
``Frame`` API.

Available tutorials
-------------------

.. table::

    +----------------------------------------------+-------------------------------------------------------------------+
    | Tutorial                                     | Description                                                       |
    +==============================================+===================================================================+
    | `Classify ice <classify-ice.rst>`_           | Load a mixed ice-water dump, print CHILL+, score cages            |
    +----------------------------------------------+-------------------------------------------------------------------+
    | `dseams <dseams.md>`_                        | CHILL+, q4/q6 and rings on NumPy arrays through ``import dseams`` |
    +----------------------------------------------+-------------------------------------------------------------------+
    | `Bulk ice <bulk-ice.rst>`_                   | CHILL+ and cages on the mixed TIP4P frame, CLI and Python         |
    +----------------------------------------------+-------------------------------------------------------------------+
    | `Ice and clathrate <ice-and-clathrate.rst>`_ | CHILL+ stars and signature cages on sI                            |
    +----------------------------------------------+-------------------------------------------------------------------+
    | `Primitive rings <rings.rst>`_               | Hydrogen-bond and cutoff ring census                              |
    +----------------------------------------------+-------------------------------------------------------------------+
    | `Order parameters <order-parameters.rst>`_   | Steinhardt q\ :sub:`l`\ and the Lechner-Dellago average           |
    +----------------------------------------------+-------------------------------------------------------------------+
    | `Confined water <confined-water.rst>`_       | In-plane RDF, monolayer rings, prisms                             |
    +----------------------------------------------+-------------------------------------------------------------------+
    | `Paper figures <paper-figures.rst>`_         | The 2020 JCIM figure calls on in-tree dumps                       |
    +----------------------------------------------+-------------------------------------------------------------------+

What a tutorial covers
----------------------

Each tutorial states prerequisites, lists learning objectives, and
walks numbered steps with expected output. The flagship path is
``ds.read`` to ``Frame.chill_plus`` to ``Frame.cages``.

Sibling packages (``seams`` in seams-core, ``require("dseams")`` in
yodaStruct) share the same engine. These tutorials stay on Python.

Getting help
------------

`How-to guides <../howto/index.rst>`_
    solve a specific task

`Python surface <../reference/python.rst>`_
    live names

`Explanation <../explanation/index.rst>`_
    why the surface is split
