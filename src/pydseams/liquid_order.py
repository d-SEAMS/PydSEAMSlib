"""Local water order from periodic neighbour distances and bond angles."""

import linkcell
import minimage
import numpy as np

def _nearest(xyz, cell, k):
    indices, squared = linkcell.knearest(xyz, cell, k)
    return np.from_dlpack(indices), np.from_dlpack(squared)


def _check_cutoff(cutoff):
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be positive and finite")


def tetrahedral_q(frame, cutoff):
    if cutoff is not None:
        _check_cutoff(cutoff)
    xyz, cell = frame._order_geometry()
    out = np.full(len(xyz), np.nan)
    if len(xyz) < 5:
        return out
    idx, squared = _nearest(xyz, cell, 4)
    valid = (idx >= 0).all(axis=1) & (squared > 0).all(axis=1)
    if cutoff is not None:
        valid &= (squared <= cutoff**2).all(axis=1)
    centers = np.flatnonzero(valid)
    if not len(centers):
        return out
    neighbors = idx[valid]
    delta = xyz[neighbors] - xyz[centers, None]
    mic = minimage.Cell.from_vectors(*cell)
    vectors = np.asarray(mic.wrap_many(delta.reshape(-1, 3))).reshape(delta.shape)
    # In a skewed cell a centered fractional image can lie outside the Wigner-Seitz cell.
    norms = np.einsum("ijk,ijk->ij", vectors, vectors)
    for row, col in zip(
        *np.where(~np.isclose(norms, squared[valid], rtol=1.0e-12, atol=1.0e-12))
    ):
        vectors[row, col] = mic.displacement_euclidean(
            xyz[centers[row]], xyz[neighbors[row, col]]
        )
    unit = vectors / np.linalg.norm(vectors, axis=2)[:, :, None]
    left, right = np.triu_indices(4, k=1)
    cosines = np.einsum("ijk,ijk->ij", unit[:, left], unit[:, right])
    out[valid] = 1.0 - 3.0 / 8.0 * np.sum((cosines + 1.0 / 3.0) ** 2, axis=1)
    return out


def local_structure_index(frame, cutoff):
    _check_cutoff(cutoff)
    xyz, cell = frame._order_geometry()
    out = np.full(len(xyz), np.nan)
    k = min(12, len(xyz) - 1)
    if k < 2:
        return out
    while True:
        _, squared = _nearest(xyz, cell, k)
        radii = np.sqrt(squared)
        inside = np.sum(radii < cutoff, axis=1)
        if np.all(inside < k) or k == len(xyz) - 1:
            break
        k = min(2 * k, len(xyz) - 1)
    gaps = np.diff(radii, axis=1)
    for count in np.unique(inside):
        if 0 < count < k:
            selected = inside == count
            out[selected] = gaps[selected, :count].var(axis=1)
    return out
