"""d-SEAMS Python front end.

``pydseams`` is the package. The compiled nanobind module is
:mod:`pydseams.yoda`. Helpers (:class:`~pydseams.frame.Frame`,
:func:`~pydseams.io.read`, ASE, solvis) sit on that surface.
``_core`` and ``cyoda`` are aliases of ``yoda``. ``Trajectory`` is an
alias of :class:`~pydseams.frame.Frame`.

Load a frame and ask for ice::

    import pydseams as ds
    frame = ds.read("water.lammpstrj")
    print(frame.chill_plus())
    print(frame.cages())

ASE Atoms work the same way::

    frame = ds.from_ase(atoms)          # default: oxygen
    atoms = frame.to_ase()
"""

from . import config
from . import yoda
from . import yoda as _core
from . import yoda as cyoda
from .frame import (
    CageScore,
    ContactPairs,
    DensityProfile,
    DomainStats,
    Frame,
    IceCounts,
    read,
)
from .io import available_readers

__version__ = "2.10.0"

# Drop-in name used in the 2.0 docs and tests
Trajectory = Frame


def from_ase(atoms, select="O", cutoff=None, bonded="auto"):
    """Build a :class:`~pydseams.frame.Frame` from an ASE ``Atoms``.

    Parameters
    ----------
    atoms : ase.Atoms
        Configuration with a nonsingular cell periodic in all three directions.
    select : str, int or sequence, optional
        Chemical symbol or atomic number of the species to analyse.
        Default ``"O"``. ``None`` keeps every atom. A sequence such as
        ``("O", "Na", "Cl")`` keeps the listed species and analyses the
        first, so ions stay in the cloud for
        :func:`pydseams.features.ion_environment`.
    cutoff : float, optional
        Neighbour cutoff in Angstroms. Default ``SEAMS_CUTOFF`` or ``3.5``.
    bonded : {"auto", "hbond", "cutoff"}, optional
        Graph for rings. ``"auto"`` uses hydrogen bonds when the
        ``Atoms`` contain H and the selected analysis cloud excludes H;
        otherwise it uses the cutoff neighbour list. Explicit ``"hbond"``
        also requires a heavy-atom selection.

    Returns
    -------
    Frame
        Analysable oxygen (or selected-species) configuration.

    Raises
    ------
    ImportError
        If ASE is not installed (``pip install 'pydseamslib[ase]'``).
    TypeError
        If ``atoms`` is not an ASE ``Atoms``.
    ValueError
        If the cell is singular, is not fully periodic, or ``select`` matches
        no atom, or if ``bonded="hbond"`` selects hydrogen.
    """
    return Frame.from_ase(
        atoms,
        select=select,
        cutoff=config.cutoff() if cutoff is None else cutoff,
        bonded=bonded,
    )


def from_arrays(positions, cell, numbers=None, cutoff=None, bonded="cutoff"):
    """Build a :class:`~pydseams.frame.Frame` from coordinates and box lengths.

    Parameters
    ----------
    positions : sequence of (x, y, z)
        Cartesian coordinates, shape ``(N, 3)``.
    cell : sequence of float
        Three orthorhombic box lengths or six LAMMPS restricted-triclinic
        bound spans and tilts.
    numbers : sequence of int, optional
        Per-atom type codes stored as ``c_type``. Default ``1`` for
        every particle.
    cutoff : float, optional
        Neighbour cutoff in Angstroms. Default ``SEAMS_CUTOFF`` or ``3.5``.
    bonded : {"auto", "hbond", "cutoff"}, optional
        Graph for rings. Default ``"cutoff"`` because this constructor
        does not attach a hydrogen cloud.

    Returns
    -------
    Frame

    Raises
    ------
    ValueError
        If ``positions`` is empty or ``cell`` is not three or six values.
    """
    return Frame.from_arrays(
        positions,
        cell,
        numbers=numbers,
        cutoff=config.cutoff() if cutoff is None else cutoff,
        bonded=bonded,
    )


def from_chemfiles(path, frame=1, **kwargs):
    """Build a :class:`~pydseams.frame.Frame` through chemfiles.

    Parameters
    ----------
    path : path-like
        Trajectory chemfiles can read (PDB, GRO, DCD, ...).
    frame : int, optional
        1-indexed frame. Default ``1``.
    **kwargs
        Forwarded to :meth:`Frame.from_chemfiles` (``cutoff``,
        ``bonded``, ``atom_type``, ``type_filter``).

    Returns
    -------
    Frame

    Raises
    ------
    RuntimeError
        If this build of seams-core did not link chemfiles.
    """
    return Frame.from_chemfiles(path, frame=frame, **kwargs)


def from_con(path, frame=1, **kwargs):
    """Build a :class:`~pydseams.frame.Frame` from an eOn ``.con`` file.

    Parameters
    ----------
    path : path-like
        eOn ``.con`` trajectory.
    frame : int, optional
        1-indexed frame. Default ``1``.
    **kwargs
        Forwarded to :meth:`Frame.from_con` (``cutoff``, ``bonded``,
        ``atom_type``).

    Returns
    -------
    Frame

    Raises
    ------
    RuntimeError
        If this build of seams-core did not link readcon-core.
    """
    return Frame.from_con(path, frame=frame, **kwargs)


def from_xyz(path, **kwargs):
    """Build a :class:`~pydseams.frame.Frame` from an XYZ file.

    Parameters
    ----------
    path : path-like
        XYZ structure.
    **kwargs
        Forwarded to :meth:`Frame.from_xyz` (``cutoff``, ``bonded``,
        ``atom_type``).

    Returns
    -------
    Frame

    Raises
    ------
    RuntimeError
        If this build of seams-core has no ``readXYZ``.
    """
    return Frame.from_xyz(path, **kwargs)


def to_solvis(frame, expand_box=True):
    """Wrap a :class:`~pydseams.frame.Frame` as a ``solvis.System``.

    Parameters
    ----------
    frame : Frame
        Configuration to view. Converted through :meth:`Frame.to_ase`.
    expand_box : bool, optional
        Passed to ``solvis.system.System``. Default ``True``.

    Returns
    -------
    solvis.system.System

    Raises
    ------
    ImportError
        If solvis is not installed (``pip install 'pydseamslib[solvis]'``).
    """
    from .solvis import to_solvis as _to

    return _to(frame, expand_box=expand_box)


__all__ = [
    "Frame",
    "Trajectory",
    "IceCounts",
    "CageScore",
    "DensityProfile",
    "ContactPairs",
    "DomainStats",
    "read",
    "from_ase",
    "from_arrays",
    "from_chemfiles",
    "from_con",
    "from_xyz",
    "to_solvis",
    "available_readers",
    "yoda",
    "_core",
    "cyoda",
    "__version__",
]
