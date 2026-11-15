project(
    'pydseamslib',
    'cpp',
    version: '2.14.0',
    default_options: [
        'warning_level=1',
        'cpp_std=c++20',
        'python.allow_limited_api=false',
    ],
    meson_version: '>=1.3.0',
)

# C++ engine. Bindings live here; seams-core is the library, not a monolith.
# libseams (the C API) and libseams_lang are shared and ship in the wheel;
# libyodaLib and the Rust crates (linkcell, minimage, readcon) are static
# inside them. readcon keeps the .con reader the 2.11 wheels had.
seams = subproject(
    'seams-core',
    default_options: [
        'default_library=static',
        'with_python=false',
        'with_tests=false',
        'with_cli=false',
        'with_linkcell=enabled',
        'with_minimage=enabled',
        'with_blas=disabled',
        'with_chemfiles=disabled',
        'with_readcon=enabled',
        'with_gpulite=disabled',
        'with_openmp_offload=disabled',
        'with_mpi=disabled',
        'with_lua=disabled',
        'with_fortran=disabled',
        'with_lammps=disabled',
        'with_gromacs=disabled',
        'with_openmm=disabled',
        'with_plumed=disabled',
        'with_ira=disabled',
        'with_spglib=disabled',
        'with_sphericart=disabled',
        'with_nauty=disabled',
        'with_csm=disabled',
    ],
)
yds_dep = seams.get_variable('yds_dep')
seams_api_lib = seams.get_variable('seams_api_lib')
seams_lang_lib = seams.get_variable('seams_lang')

# seams.h, and the seams_version.h it includes from the subproject build
# tree. seams-core declares no dependency object for libseams, and
# include_directories() into a subproject is a sandbox violation.
dseams_api_dep = declare_dependency(
    link_with: seams_api_lib,
    compile_args: [
        '-I' + meson.project_source_root() / 'subprojects' / 'seams-core' / 'include',
        '-I' + meson.project_build_root() / 'subprojects' / 'seams-core' / 'src',
    ],
)

nanobind_dep = dependency(
    'nanobind',
    default_options: [
        'python.allow_limited_api=' + get_option('python.allow_limited_api').to_string(),
    ],
)

py = import('python').find_installation(pure: false)
# Stable-ABI wheel builds opt in; source builds use the interpreter ABI.
limited_api = py.language_version().version_compare('>=3.12') ? '3.12' : ''

_link_args = []
if host_machine.system() == 'linux'
    _link_args += ['-Wl,--enable-new-dtags']
endif

py.extension_module(
    'yoda',
    'src/bindings.cpp',
    dependencies: [yds_dep, nanobind_dep],
    link_args: _link_args,
    install: true,
    subdir: 'pydseams',
    limited_api: limited_api,
)

py.install_sources(
    [
        'src/pydseams/__init__.py',
        'src/pydseams/band.py',
        'src/pydseams/features.py',
        'src/pydseams/frame.py',
        'src/pydseams/aseio.py',
        'src/pydseams/io.py',
        'src/pydseams/solvis.py',
        'src/pydseams/config.py',
    ],
    pure: false,
    subdir: 'pydseams',
)

py.install_sources(
    [
        'src/pydseams/adapters/__init__.py',
        'src/pydseams/adapters/_ase.py',
        'src/pydseams/adapters/_mdanalysis.py',
        'src/pydseams/adapters/_ovito.py',
    ],
    pure: false,
    subdir: 'pydseams/adapters',
)

py.install_sources(
    ['src/pydseamslib/__init__.py'],
    pure: false,
    subdir: 'pydseamslib',
)

# The 2.14 API: import dseams.
subdir('src/dseams')

# The sdist keeps the patched linkcell and the vesin and gpulite overlays.
meson.add_dist_script(py, files('scripts/dist_subprojects.py'))

# Development-tree tests: `meson test -C <builddir>`. The wheel and CI
# run the same files with pytest against the installed package.
_engine_src = meson.project_source_root() / 'subprojects' / 'seams-core'
_dseams_env = environment()
_dseams_env.set('DSEAMS_BUILD', meson.project_build_root())
_dseams_env.set('DSEAMS_SRC', meson.project_source_root())
_dseams_env.set('SEAMS_CORE_SRC', _engine_src)
_dseams_env.set('SEAMS_LIB', seams_api_lib.full_path())
_dseams_env.set('SEAMS_LANG_LIB', seams_lang_lib.full_path())
test(
    'python_dseams',
    py,
    args: [files('tests/dseams/run_unittest.py')],
    env: _dseams_env,
    depends: [dseams_core, seams_lang_lib],
    timeout: 180,
)
_ctypes_env = environment()
_ctypes_env.set('PYTHONPATH', meson.project_source_root() / 'src')
_ctypes_env.set('SEAMS_LIB', seams_api_lib.full_path())
_ctypes_env.set('SEAMS_CORE_SRC', _engine_src)
test(
    'python_seams_ctypes',
    py,
    args: [
        '-m',
        'pytest',
        '-q',
        '-p',
        'no:cacheprovider',
        '--noconftest',
        meson.project_source_root() / 'tests' / 'dseams' / 'test_seams_package.py',
        '--doctest-modules',
        meson.project_source_root() / 'src' / 'seams',
    ],
    env: _ctypes_env,
    depends: [seams_api_lib],
    timeout: 180,
)
_lang_env = environment()
_lang_env.set('SEAMS_LANG_LIB', seams_lang_lib.full_path())
test(
    'python_chill_plus_ctypes',
    py,
    args: [files('tests/dseams/test_chill_plus_ctypes.py')],
    env: _lang_env,
    depends: [seams_lang_lib],
    timeout: 60,
)
_examples_env = environment()
_examples_env.set('SEAMS_CORE_SRC', _engine_src)
_examples_env.set(
    'SEAMS_LIBDIR',
    meson.project_build_root() / 'subprojects' / 'seams-core' / 'src',
)
test(
    'docs_api_python',
    py,
    args: [
        files('examples/python/python_api.py'),
        meson.project_build_root() / 'subprojects' / 'seams-core' / 'src',
    ],
    env: _examples_env,
    depends: [seams_api_lib],
    timeout: 60,
)
test(
    'docs_minimal_python',
    py,
    args: [
        files('examples/python/minimal_py.py'),
        _engine_src / 'input' / 'traj' / 'genice_sI.lammpstrj',
        meson.project_build_root() / 'subprojects' / 'seams-core' / 'src',
    ],
    env: _examples_env,
    depends: [seams_api_lib],
    timeout: 60,
)
test(
    'docs_core_counts',
    py,
    args: [files('examples/python/core_counts.py')],
    env: _examples_env,
    depends: [seams_api_lib],
    timeout: 120,
)
