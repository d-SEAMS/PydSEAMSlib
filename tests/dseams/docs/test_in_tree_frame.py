# Copyright (c) d-SEAMS developers
"""The tutorial job imports the 2.14 frame.

PyPI ``pydseamslib`` 2.11.1 has no ``dseams.userframe``, its ``chill_plus``
has no ``graph`` argument, and it labels the mixed TIP4P dump 12
interClathrate and 238 water. This file fails when that package is the
one on ``sys.path``.
"""

import inspect


def test_frame_is_the_in_tree_214_package():
    import dseams
    import pydseams
    from dseams import userframe

    assert dseams.__version__ == "2.14.1"
    assert pydseams.__version__ == "2.14.1"
    assert "graph" in inspect.signature(userframe.Frame.chill_plus).parameters
    assert "edge" in inspect.signature(userframe.Frame.cages_by_signature).parameters
