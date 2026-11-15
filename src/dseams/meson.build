# Nanobind module over libseams, plus the pure-Python package.
# The top-level meson.build sets py, nanobind_dep, limited_api and
# dseams_api_dep (libseams from the seams-core subproject).

dseams_core = py.extension_module(
  '_core',
  files('_core.cpp'),
  dependencies: [nanobind_dep, dseams_api_dep],
  link_args: _link_args,
  install: true,
  subdir: 'dseams',
  limited_api: limited_api)
py.install_sources(
  files(
    '__init__.py',
    '__main__.py',
    'frames.py',
    'userframe.py',
  ),
  subdir: 'dseams',
  pure: false)
