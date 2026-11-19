======================
View a frame in solvis
======================



Problem
-------

You have a classified ``Frame`` and want a ``solvis.System`` for PyVista
visualisation.

Install the extra
-----------------

.. code:: bash

    pip install 'pydseamslib[solvis]'

This extra installs ASE and ``solvis-tools``. Without it, ``to_solvis``
raises ``ImportError``.

Wrap the frame
--------------

Classify the frame before exporting so ``ice_type``, ``hc``, and ``ddc``
arrays are present on the wrapped ``Atoms``:

.. code:: python

    import pydseams as ds

    a = 2.75 * 4 / 3**0.5
    basis = (
        (0, 0, 0),
        (0.5, 0.5, 0),
        (0.5, 0, 0.5),
        (0, 0.5, 0.5),
        (0.25, 0.25, 0.25),
        (0.75, 0.75, 0.25),
        (0.75, 0.25, 0.75),
        (0.25, 0.75, 0.75),
    )
    xyz = [
        [(i + x) * a, (j + y) * a, (k + z) * a]
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for x, y, z in basis
    ]
    frame = ds.from_arrays(xyz, [3 * a, 3 * a, 3 * a])
    frame.chill_plus()
    frame.cages()

    system = frame.to_solvis()
    # Equivalent: system = ds.to_solvis(frame)
    assert type(system).__name__ == "System"

``solvis`` wraps the ``Atoms`` object returned by ``Frame.to_ase()``, which
carries ``arrays["ice_type"]`` after CHILL/CHILL+ and ``arrays["hc"]`` /
``arrays["ddc"]`` after ``cages()``. Pass ``expand_box=False`` to disable
solvis box expansion (default ``True``).

Interactive plots and screenshots follow the
`solvis examples <https://github.com/amritagos/solvis>`_ (press ``q`` to
close a PyVista window).

See also
--------

`Classify ASE Atoms <ase.rst>`_
    ``from_ase`` / ``to_ase``

`The yoda surface <../explanation/yoda-surface.rst>`_
    why visualisation stays in Python
