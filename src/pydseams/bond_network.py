"""Connected components and residence statistics for sampled bond networks."""

from operator import index

import numpy as np


def components(n_sites, pairs):
    """Return per-site component labels and component sizes.

    Each bond joins two zero-based site indices regardless of direction.
    Labels follow the first site encountered in each component. Isolated
    sites have size one.
    """
    n_sites = index(n_sites)
    if n_sites < 0:
        raise ValueError("n_sites must be nonnegative")
    parent = list(range(n_sites))
    size = [1] * n_sites

    def root(site):
        while parent[site] != site:
            parent[site] = parent[parent[site]]
            site = parent[site]
        return site

    for left, right in pairs:
        left, right = index(left), index(right)
        if not (0 <= left < n_sites and 0 <= right < n_sites):
            raise ValueError("bond site index is outside the network")
        a, b = root(left), root(right)
        if a != b:
            if size[a] < size[b]:
                a, b = b, a
            parent[b] = a
            size[a] += size[b]

    labels_by_root = {}
    labels = np.empty(n_sites, dtype=np.intp)
    for site in range(n_sites):
        representative = root(site)
        labels[site] = labels_by_root.setdefault(representative, len(labels_by_root))
    sizes = np.fromiter((size[r] for r in labels_by_root), dtype=np.intp)
    return labels, sizes


def residence_counts(membership, lags):
    """Return lag, at-risk origin, continuous and intermittent counts.

    The boolean matrix has rows at equally spaced saved frames and columns
    for stable site identities. An origin counts only when its endpoint
    exists. Continuous residence requires membership in every saved frame
    in the interval. Intermittent residence requires both endpoints.
    """
    membership = np.asarray(membership)
    if membership.ndim != 2 or membership.dtype != np.dtype(bool):
        raise ValueError("membership must be a two-dimensional boolean array")
    n_frames, n_sites = membership.shape
    lags = [index(lag) for lag in lags]
    if any(not 0 <= lag < n_frames for lag in lags):
        raise ValueError("each lag must be within the saved frames")

    streaks = np.empty(membership.shape, dtype=np.intp)
    running = np.zeros(n_sites, dtype=np.intp)
    for frame in range(n_frames - 1, -1, -1):
        running = np.where(membership[frame], running + 1, 0)
        streaks[frame] = running

    counts = np.empty((len(lags), 4), dtype=np.int64)
    for row, lag in enumerate(lags):
        stop = n_frames - lag
        start = membership[:stop]
        counts[row] = (
            lag,
            np.count_nonzero(start),
            np.count_nonzero(streaks[:stop] > lag),
            np.count_nonzero(start & membership[lag:]),
        )
    return counts
