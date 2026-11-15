Python
======

The frame is a context manager. The ring histogram and the CHILL+ labels
are NumPy arrays the library fills in place.

.. literalinclude:: ../../../examples/python/lang_tasks_python.py
   :language: python

``meson test docs_lang_python`` runs this file, and
``tests/dseams/test_doc_snippets.py`` checks its output on
``input/traj/ice1h.lammpstrj`` against the installed wheel.
