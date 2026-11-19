=========
Tutorials
=========


Tutorials are learning-oriented. Follow them in order to learn the
``Frame`` and ``dseams`` Python APIs.

Available tutorials
-------------------

.. table::

    +----------------------------------------------+--------------------------------------------------------------------------------+
    | Tutorial                                     | Description                                                                    |
    +==============================================+================================================================================+
    | `Classify ice <classify-ice.rst>`_           | Load a mixed ice-water dump, print CHILL+, score cages                         |
    +----------------------------------------------+--------------------------------------------------------------------------------+
    | `dseams <dseams.md>`_                        | CHILL+, q4/q6, rings, and cages on NumPy arrays through ``import dseams``      |
    +----------------------------------------------+--------------------------------------------------------------------------------+
    | `Bulk ice <bulk-ice.rst>`_                   | Compare cutoff CHILL+, mutual-4 CHILL+, and cages on mixed TIP4P and cubic ice |
    +----------------------------------------------+--------------------------------------------------------------------------------+
    | `Ice and clathrate <ice-and-clathrate.rst>`_ | CHILL+ stars and signature cages on sI clathrate                               |
    +----------------------------------------------+--------------------------------------------------------------------------------+
    | `Primitive rings <rings.rst>`_               | Hydrogen-bond and cutoff primitive ring census                                 |
    +----------------------------------------------+--------------------------------------------------------------------------------+
    | `Order parameters <order-parameters.rst>`_   | Steinhardt q\ :sub:`l`\ and the Lechner-Dellago average                        |
    +----------------------------------------------+--------------------------------------------------------------------------------+
    | `Confined water <confined-water.rst>`_       | In-plane RDF, monolayer rings, and prisms                                      |
    +----------------------------------------------+--------------------------------------------------------------------------------+
    | `Paper figures <paper-figures.rst>`_         | Reproduce key d-SEAMS paper structural analyses in Python                      |
    +----------------------------------------------+--------------------------------------------------------------------------------+

What a tutorial covers
----------------------

Each tutorial states prerequisites, lists learning goals, and walks
through runnable steps with expected output. The main workflow is
``ds.read`` (or ``ds.from_arrays``) to ``Frame.chill_plus`` and
``Frame.cages``.

Getting help
------------

`How-to guides <../howto/index.rst>`_
    solve a specific task

`Python surface <../reference/python.rst>`_
    live names and records

`Explanation <../explanation/index.rst>`_
    why the surface is split
