C-ABI Ice Ih census
===================

For environments calling ``libseams_lang`` directly over ``ctypes``
without the compiled Python extensions, :class:`Frame` below wraps the C
frame handle as a context manager and fills NumPy arrays in place.
Standard Python workflows installed via ``pip install pydseamslib`` can
call :mod:`dseams` or :class:`dseams.userframe.Frame` directly instead.

.. literalinclude:: ../../../examples/python/lang_tasks_python.py
   :language: python
