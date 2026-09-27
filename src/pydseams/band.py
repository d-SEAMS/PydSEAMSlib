"""Pairs that enter or leave the cutoff graph along a sequence of frames.

A primitive-ring count labels a frame. It orders a path when it moves
one way. ``band_rows`` names the pairs that enter or leave the graph
and records the shortest distance inside each atom type.
"""

from __future__ import annotations

from . import yoda


def bond_pairs(frame):
    """Undirected pairs in the cutoff graph, as index tuples ``(i, j)`` with ``i < j``."""
    pairs = set()
    for i, row in enumerate(frame.bonds_by_index):
        for j in row[1:]:
            if j > i:
                pairs.add((i, j))
    return pairs


def _distance(frame, i, j):
    a = frame.cloud.pts[i]
    b = frame.cloud.pts[j]
    return ((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2) ** 0.5


def shortest_within_type(frame, type_code):
    """Shortest free-space distance between two atoms of ``type_code``.

    Returns ``None`` when fewer than two atoms carry that type.
    """
    indices = [
        i for i, pt in enumerate(frame.cloud.pts) if int(pt.c_type) == int(type_code)
    ]
    best = None
    for a in range(len(indices)):
        for b in range(a + 1, len(indices)):
            distance = _distance(frame, indices[a], indices[b])
            if best is None or distance < best:
                best = distance
    return best


def ring_census(frame, max_size=8):
    """Count primitive rings by size, up to ``max_size``."""
    counts = {}
    for ring in yoda.ringNetwork(frame.bonds_by_index, int(max_size)):
        size = len(ring)
        counts[size] = counts.get(size, 0) + 1
    return counts


def band_rows(frames, max_ring_size=8):
    """One row per frame: bonds that entered or left, shortest same-type distances, rings.

    Parameters
    ----------
    frames : sequence of Frame
        Images in path order. Each frame uses its own cutoff graph.
    max_ring_size : int, optional
        Largest primitive ring. Default 8.

    Returns
    -------
    list of dict
        ``index``, ``n_bonds``, ``entered``, ``left``,
        ``shortest_within_type``, ``rings``. The first row has empty
        ``entered`` and ``left``.
    """
    rows = []
    previous = None
    for index, frame in enumerate(frames):
        pairs = bond_pairs(frame)
        if previous is None:
            entered, left = [], []
        else:
            entered = sorted(pairs - previous)
            left = sorted(previous - pairs)
        types = sorted({int(pt.c_type) for pt in frame.cloud.pts})
        rows.append(
            {
                "index": index,
                "n_bonds": len(pairs),
                "entered": entered,
                "left": left,
                "shortest_within_type": {
                    t: shortest_within_type(frame, t) for t in types
                },
                "rings": ring_census(frame, max_ring_size),
            }
        )
        previous = pairs
    return rows
