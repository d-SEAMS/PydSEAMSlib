=========
Why Frame
=========


CHILL, rings, cages, Steinhardt parameters, and SOAP descriptors are
all per-snapshot quantities. ``Frame`` models one atomic configuration:
its coordinates, simulation box, neighbour lists, bonded graph, and
per-atom classifications. ``Frame.cages`` returns the topological cage
score (hexagonal cages = ice Ih, double-diamond cages = ice Ic,
neither = water).

Cutoff neighbour lists (``neighbor_list``, CHILL, cutoff-graph cages)
use vesin when linked, whereas seeded hysteresis cages
(``cages(seeded=True)``) use linkcell k-nearest neighbours. Those
neighbour backends serve distinct graph definitions and are not
interchangeable.

Runtime configuration
---------------------

``pydseams.config`` shares the twelve-factor configuration precedence
used across d-SEAMS: module defaults are overridden by ``SEAMS_CONFIG``
or ``./seams.env``, which are overridden by process environment variables,
which are overridden by explicit function arguments. ``Frame`` reads
``SEAMS_FRAME`` and ``SEAMS_CUTOFF`` when ``frame`` or ``cutoff`` are omitted.

One configuration at a time
---------------------------

``Frame`` holds a single snapshot. Multi-frame workflows step through a
trajectory explicitly with ``read(path, frame=n)`` (1-indexed) or
``frame.load_frame(n)`` on a file-backed ``Frame``. Frames constructed
from in-memory arrays or ASE ``Atoms`` have no backing file path and do
not support ``load_frame``.

``Trajectory`` is an alias
--------------------------

``Trajectory`` is a compatibility alias of ``Frame``. It holds one active
frame at a time and loads subsequent frames via ``load_frame(n)`` when
iterated (for example inside ``IceFeaturizer.transform``).
