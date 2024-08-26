project(
    'PydSEAMSlib',
    'cpp',
    version: '1.0.1',
    default_options: ['warning_level=1', 'cpp_std=c++14'],
)

host_system = host_machine.system()

# Add C++ compiler options
# Note that fast-math gives incorrect results, and the third set is for debugging
_args = ['-funroll-loops']  # ['-ffast-math', '-funroll-loops'] #['-ggdb', '-Og', '-fno-inline-functions']
_deps = []  # Dependencies
_linkto = []  # All the sub-libraries
_incdirs = []  # All the includes

# ---------------------- Library Dependencies
ydslib = subproject(
    'seams-core',
    default_options: ['default_library=static', 'wrap_mode=forcefallback'],
).get_variable(
    'yds_dep',
)
_deps += ydslib

# ---------------------- Bindings
py_mod = import('python')
py = py_mod.find_installation(pure: false)
python_dep = py.dependency()
pyb11f_deps = [python_dep, dependency('pybind11')]

_deps += [pyb11f_deps, declare_dependency(link_args: '-lstdc++')]

# cyoda, for the bindings to seams-core
py.extension_module(
    'cyoda',
    sources: ['pyb11_srcs/py_one.cc'],
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
    subdir: 'pydseamslib/',
)

# pydseamslib main package
py.install_sources(
    ['src/__init__.py'],
    pure: false,  # install next to compiled extension
    subdir: 'pydseamslib/',
)

py.install_sources(
    ['src/adapters/__init__.py', 'src/adapters/_ase.py'],
    pure: false,
    subdir: 'pydseamslib/adapters',
)
