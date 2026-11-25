=======================================
Classify ice on the mixed TIP4P example
=======================================


You classify a single frame two ways: CHILL+ (a four-neighbour
staggered/eclipsed star per oxygen) and topological cages (HC/DDC
membership). Comparing both methods on a mixed frame and a cubic
crystal shows why the two scores can disagree.

Prerequisites
-------------

- Python 3.11+ with ``pydseamslib`` installed (``pip install 'pydseamslib==2.14.2'``)

- Sample trajectories ``exampleTraj.lammpstrj`` (250 TIP4P waters;
  LAMMPS type 2 oxygen, type 1 hydrogen) and ``mW_cubic.lammpstrj``
  (4096 coarse-grained mW sites, LAMMPS type 1)

.. code:: bash

    mkdir -p input/traj
    curl -sLo input/traj/exampleTraj.lammpstrj \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/v2.14.2/input/traj/exampleTraj.lammpstrj
    curl -fsSL -o input/fetch.py \
      https://raw.githubusercontent.com/d-SEAMS/seams-core/v2.14.2/input/fetch.py
    python3 input/fetch.py mW_cubic.lammpstrj

Learning goals
--------------

In this tutorial you learn to:

1. Read a LAMMPS dump with ``ds.read`` and select the oxygen type.

2. Run CHILL+ on both the cutoff neighbour list and the mutual
   four-nearest neighbour list.

3. Run topological cages on the same oxygens and inspect HC / DDC
   membership.

4. Explain why CHILL+ and cage scores differ on a mixed frame yet agree
   on a cubic crystal.

Step 1: Load the mixed TIP4P frame and run cutoff CHILL+
--------------------------------------------------------

CHILL+ inspects each oxygen's bonded neighbours. The default
``Frame.chill_plus()`` builds a cutoff neighbour list (``cutoff=3.5`` A),
classifying each bond as staggered or eclipsed from the pair
correlation :math:`c_{ij}`:

.. table::

    +----------------+---------------------------------------------+
    | label          | four-neighbour star                         |
    +================+=============================================+
    | cubic          | four staggered                              |
    +----------------+---------------------------------------------+
    | hexagonal      | three staggered and one eclipsed            |
    +----------------+---------------------------------------------+
    | clathrate      | four eclipsed                               |
    +----------------+---------------------------------------------+
    | interClathrate | three eclipsed                              |
    +----------------+---------------------------------------------+
    | interfacial    | ice-like, not one of the above              |
    +----------------+---------------------------------------------+
    | water          | fewer than four bonds, or none of the above |
    +----------------+---------------------------------------------+

.. code:: python

    import pydseams as ds

    frame = ds.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
    print("n_atoms:", frame.n_atoms, "atom_type:", frame.atom_type)

    counts_cutoff = frame.chill_plus()
    print("cutoff CHILL+:", counts_cutoff)
    assert str(counts_cutoff) == "IceCounts(interClathrate=9, water=241)"

Expected output:

::

    n_atoms: 250 atom_type: 2
    cutoff CHILL+: IceCounts(interClathrate=9, water=241)

Step 2: Compare with the mutual four-nearest CHILL+ graph
---------------------------------------------------------

``dseams.userframe.Frame.chill_plus`` accepts a ``graph`` argument. With
``graph="mutual"``, it keeps only mutual four-nearest neighbours
(``k = 4``) instead of every pair inside the 3.5 A cutoff:

.. code:: python

    from dseams import userframe

    uframe = userframe.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
    counts_mutual = uframe.chill_plus(graph="mutual")
    print("mutual CHILL+:", counts_mutual)
    assert str(counts_mutual) == "IceCounts(clathrate=1, interClathrate=7, water=242)"

Expected output:

::

    mutual CHILL+: IceCounts(clathrate=1, interClathrate=7, water=242)

.. table::

    +----------------+--------------------+
    | label          | atoms (mutual k=4) |
    +================+====================+
    | clathrate      |                  1 |
    +----------------+--------------------+
    | interClathrate |                  7 |
    +----------------+--------------------+
    | water          |                242 |
    +----------------+--------------------+

One oxygen has four eclipsed bonds and seven have three eclipsed bonds;
the rest are water. Neither star forms a cubic or hexagonal ice motif.

Step 3: Score cages on the same frame
-------------------------------------

Cages ask a topological question: does this oxygen belong to a complete
hexagonal cage (HC, ice Ih, 12 waters) or a complete double-diamond
cage (DDC, ice Ic, 14 waters)? The search walks primitive six-rings on
the bonded graph (default seeded hysteresis):

.. code:: python

    cages = frame.cages()
    print("cages:", cages)
    assert str(cages) == "CageScore(n_ih=0, n_ic=0, n_water=250)"

Expected output:

::

    cages: CageScore(n_ih=0, n_ic=0, n_water=250)

This mixed frame has no complete HC or DDC.

Step 4: Why the two scores differ
---------------------------------

CHILL+ is a local four-neighbour star. A molecule is labelled ice-like
whenever its immediate four bonds match a tetrahedral motif, even if
those neighbours do not close into a polyhedral cage.

Cages implement topological unit matching (TUM). A molecule is ice Ih
only when it participates in a closed hexagonal cage, and ice Ic only
when it participates in a closed double-diamond cage. Partial local
order without cage closure remains water.

Step 5: Contrast with a cubic crystal
-------------------------------------

On the 4096-molecule cubic mW crystal (``mW_cubic.lammpstrj``, type 1),
every molecule has four staggered bonds and sits inside a complete DDC,
so CHILL+ and cage membership both classify all 4096 atoms as cubic:

.. code:: python

    cubic = ds.read("input/traj/mW_cubic.lammpstrj", frame=1, atom_type=1)
    cubic_chill = cubic.chill_plus()
    cubic_cages = cubic.cages()
    print("cubic CHILL+:", cubic_chill)
    print("cubic cages:", cubic_cages)
    assert cubic_chill.cubic == 4096
    assert cubic_cages.n_ic == 4096

Expected output:

::

    cubic CHILL+: IceCounts(cubic=4096)
    cubic cages: CageScore(n_ih=0, n_ic=4096, n_water=0)

Troubleshooting
---------------

- Empty histogram or ``ValueError`` on types: verify ``atom_type``. Mixed
  TIP4P uses type 2 for oxygen; single-site mW uses type 1.

- ``TypeError`` passing ``graph=`` to ``pydseams.Frame.chill_plus``: use
  ``dseams.userframe.Frame.chill_plus(graph="mutual")`` (or
  ``dseams.chill_plus(..., graph=dseams.MUTUAL)``) when switching
  between cutoff and mutual-4 graphs.

Next steps
----------

`Ice and clathrate <ice-and-clathrate.rst>`_
    clathrate cage signatures

`Primitive rings <rings.rst>`_
    ring census on hydrogen-bond vs cutoff graphs

`dseams <dseams.md>`_
    NumPy array interface via ``import dseams``
