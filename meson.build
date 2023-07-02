project('pyseams', 'cpp',
       version: '1.0.1',
    default_options: ['warning_level=3', 'cpp_std=c++17'])

host_system = host_machine.system()

# Add C++ compiler options
# Note that fast-math gives incorrect results, and the third set is for debugging
   _args = ['-funroll-loops'] # ['-ffast-math', '-funroll-loops'] #['-ggdb', '-Og', '-fno-inline-functions']
   _deps = [] # Dependencies
   _linkto = [] # All the sub-libraries
   _incdirs = [] # All the includes

# ---------------------- Library Dependencies
seams_core_deps = subproject('seams-core', version: '1.0.1').get_variable('_deps')
seams_core_inc = subproject('seams-core', version: '1.0.1').get_variable('_incdirs')
seams_core_linkto = subproject('seams-core', version: '1.0.1').get_variable('_linkto')
ydslib = subproject('seams-core', version: '1.0.1').get_variable('ydslib')
seams_core_dep = declare_dependency(link_with: [seams_core_linkto, ydslib],
                                    dependencies: seams_core_deps,
                                    include_directories: seams_core_inc)
_deps += seams_core_dep

# ---------------------- Bindings
py_mod = import('python')
py = py_mod.find_installation(pure: false)
python_dep = py.dependency()
pyb11f_deps = [ python_dep, dependency('pybind11') ]

_deps += [pyb11f_deps,
          declare_dependency(link_args: '-lstdc++') ]

# cyoda, for the bindings to seams-core
py.extension_module(
  'cyoda',
  sources : [
    'pyb11_srcs/py_one.cc',
  ],
  dependencies: _deps,
  link_with: _linkto,
  cpp_args : _args,
  install: true,
)

