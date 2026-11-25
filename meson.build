project(
    'pydseamslib',
    'cpp',
    # [[[cog
    # import tomllib
    # version = tomllib.load(open("pyproject.toml", "rb"))["project"]["version"]
    # cog.outl(f"    version: '{version}',")
    # ]]]
    version: '2.14.2',
    # [[[end]]]
    default_options: [
        'warning_level=1',
        'cpp_std=c++20',
        'python.allow_limited_api=false',
    ],
    meson_version: '>=1.3.0',
)

# Readers. auto is seams-core's own detection: chemfiles from pkg-config
# or the patched library that seams-core builds, readcon-core from
# pkg-config or the readcon-core wrap. enabled fails the build when the
# library is missing. disabled keeps LAMMPS and XYZ only. When a reader
# is linked, seams-core puts -DSEAMS_HAS_CHEMFILES / -DSEAMS_HAS_READCON
# on yds_dep, and this extension compiles against that dependency.
chemfiles_choice = 'auto'
if get_option('with_chemfiles').enabled()
    chemfiles_choice = 'enabled'
elif get_option('with_chemfiles').disabled()
    chemfiles_choice = 'disabled'
endif
readcon_choice = 'auto'
if get_option('with_readcon').enabled()
    readcon_choice = 'enabled'
elif get_option('with_readcon').disabled()
    readcon_choice = 'disabled'
endif

# C++ engine. An installed libseams (pkg-config seams, seams-core and
# seams-lang) is used when it is on the prefix. Otherwise seams-core is
# the subproject, and libyodaLib plus the Rust crates are static inside
# the shared libraries that the wheel ships.
seams_core_dep = dependency('seams-core', version: '>=2.14.2', required: false)
seams_dep = dependency('seams', version: '>=2.14.2', required: false)
seams_lang_dep = dependency('seams-lang', version: '>=2.14.2', required: false)
use_installed_engine = (
    seams_core_dep.found()
    and seams_dep.found()
    and seams_lang_dep.found()
)
if use_installed_engine
    yds_dep = seams_core_dep
    dseams_api_dep = seams_dep
    seams_api_lib = seams_dep
    seams_lang_lib = seams_lang_dep
else
    seams = subproject(
        'seams-core',
        default_options: [
            'default_library=static',
            'with_tests=false',
            'with_cli=false',
            'with_linkcell=enabled',
            'with_minimage=enabled',
            'with_blas=disabled',
            'with_chemfiles=' + chemfiles_choice,
            'with_readcon=' + readcon_choice,
            'with_gpulite=disabled',
            'with_openmp_offload=disabled',
            'with_mpi=disabled',
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
    # seams.h, and the seams_version.h it includes from the subproject
    # build tree. include_directories() into a subproject is a sandbox
    # violation, so the include path is spelled out.
    dseams_api_dep = declare_dependency(
        link_with: seams_api_lib,
        compile_args: [
            '-I' + meson.project_source_root() / 'subprojects' / 'seams-core' / 'include',
            '-I' + meson.project_build_root() / 'subprojects' / 'seams-core' / 'src',
        ],
    )
endif

# seams-core copies the vendored chemfiles library into libdir with no
# install tag. A wheel built with --tags=runtime,python-runtime then
# omits it. Install the file the loader opens, named by patchelf, otool,
# or objdump, with the runtime tag. auditwheel and delocate repair that
# copy. An installed engine is found by the loader and is not staged here.
# with_chemfiles=enabled, or a subproject that actually vendored the
# library, fails configure when the tool or the library is missing.
if not use_installed_engine and not get_option('with_chemfiles').disabled()
    _chem_libdir = meson.project_build_root() / 'subprojects' / 'seams-core' / 'chemfiles-prefix' / 'lib'
    _chem_vendor = import('fs').is_dir(_chem_libdir)
    if get_option('with_chemfiles').enabled() or _chem_vendor
        if host_machine.system() == 'darwin'
            _soname_tool = find_program(
                'otool',
                required: get_option('with_chemfiles'),
            )
        else
            _soname_tool = find_program('patchelf', required: false)
            if not _soname_tool.found()
                _soname_tool = find_program(
                    'objdump',
                    required: get_option('with_chemfiles'),
                )
            endif
        endif
        if get_option('with_chemfiles').enabled() and not _chem_vendor
            error(
                'with_chemfiles is enabled but the vendored chemfiles library was not built',
            )
        endif
        if not _soname_tool.found()
            if host_machine.system() == 'darwin'
                error('reading the vendored chemfiles install name needs otool')
            else
                error(
                    'reading the vendored chemfiles SONAME needs patchelf or objdump',
                )
            endif
        endif
        _chem_probe = run_command(
            find_program('sh'),
            '-c',
            '''
libdir=$1
tool=$2
real=
best=0
for f in "$libdir"/libseams-chemfiles*; do
  [ -f "$f" ] && [ ! -L "$f" ] || continue
  case "$f" in
    *.a) continue ;;
  esac
  bytes=$(wc -c < "$f")
  if [ "$bytes" -gt "$best" ]; then
    best=$bytes
    real=$f
  fi
done
if [ -z "$real" ]; then
  echo "vendored chemfiles library not found in $libdir" >&2
  exit 1
fi
base=${tool##*/}
if [ "$base" = patchelf ]; then
  name=$("$tool" --print-soname "$real")
elif [ "$base" = otool ]; then
  name=$("$tool" -D "$real" | awk 'NR==2 { print $1 }')
  name=${name##*/}
else
  name=$("$tool" -p "$real" | awk '/SONAME/ { print $2; exit }')
fi
if [ -z "$name" ]; then
  echo "could not read the SONAME of $real with $tool" >&2
  exit 1
fi
printf '%s\n%s\n' "$real" "$name"
''',
            'chemfiles-soname',
            _chem_libdir,
            _soname_tool.full_path(),
            check: true,
        )
        _chem_lines = _chem_probe.stdout().strip().split('\n')
        if _chem_lines.length() != 2 or _chem_lines[0] == '' or _chem_lines[1] == ''
            error('vendored chemfiles SONAME probe returned no library name')
        endif
        custom_target(
            'bundle-seams-chemfiles',
            output: _chem_lines[1],
            command: ['cp', _chem_lines[0], '@OUTPUT@'],
            install: true,
            install_dir: get_option('libdir'),
            install_tag: 'runtime',
            install_mode: 'rwxr-xr-x',
            build_by_default: true,
        )
    endif
endif

py = import('python').find_installation(pure: false)
# An installed nanobind module ships headers and src/nb_combined.cpp, not a
# library meson can find. Copy those sources into the build directory and
# compile them. A wheel build with no module still uses the wrap.
nb_root_r = run_command(
    py,
    '-c',
    'import nanobind, pathlib; print(pathlib.Path(nanobind.__file__).resolve().parent)',
    check: false,
)
nb_names = []
if nb_root_r.returncode() == 0 and nb_root_r.stdout().strip() != ''
    nb_list_r = run_command(
        py,
        '-c',
        'import pathlib, sys; root = pathlib.Path(sys.argv[1]); print(" ".join(sorted(p.name for p in root.iterdir() if p.is_file())))',
        nb_root_r.stdout().strip() / 'src',
        check: false,
    )
    if nb_list_r.returncode() == 0 and nb_list_r.stdout().strip() != ''
        nb_names = nb_list_r.stdout().split()
    endif
endif
if nb_names.length() > 0 and nb_names.contains('nb_combined.cpp')
    nb_root = nb_root_r.stdout().strip()
    nb_inc = nb_root / 'include'
    nb_robin = nb_root / 'ext' / 'robin_map' / 'include'
    nb_rest = []
    foreach nb_name : nb_names
        if nb_name != 'nb_combined.cpp'
            nb_rest += nb_name
        endif
    endforeach
    # First output is the translation unit. The other files sit beside it
    # so its quoted includes resolve, and they are not compiled on their own.
    nb_outputs = ['nb_combined.cpp'] + nb_rest
    nb_staged = custom_target(
        'stage-nanobind',
        output: nb_outputs,
        command: [
            py,
            files('scripts/stage_nanobind.py'),
            nb_root / 'src',
            '@OUTDIR@',
            nb_outputs,
        ],
    )
    nb_compile = ['-I' + nb_inc, '-I' + nb_robin]
    nb_link = []
    if meson.get_compiler('cpp').get_argument_syntax() != 'msvc'
        if host_machine.system() == 'darwin'
            nb_compile += ['-fno-strict-aliasing']
        else
            nb_compile += [
                '-fno-strict-aliasing',
                '-ffunction-sections',
                '-fdata-sections',
            ]
            nb_link += ['-Wl,--gc-sections']
        endif
    endif
    if py.get_variable('Py_GIL_DISABLED', 0) != 0
        nb_compile += ['-DNB_FREE_THREADED']
    endif
    if not get_option('debug')
        nb_compile += ['-DNB_COMPACT_ASSERTIONS']
        if host_machine.system() == 'darwin'
            nb_link += ['-Wl,-dead_strip', '-Wl,-x', '-Wl,-S']
        elif host_machine.system() != 'windows'
            nb_link += ['-Wl,-s']
        endif
    endif
    if host_machine.system() == 'darwin'
        py_soabi = py.get_variable('SOABI')
        if py_soabi.startswith('cp')
            resp_file = nb_root / 'cmake' / 'darwin-ld-cpython.sym'
        else
            resp_file = nb_root / 'cmake' / 'darwin-ld-pypy.sym'
        endif
        nb_link += ['-Wl,@' + resp_file]
    endif
    limited_api_args = []
    if get_option('python.allow_limited_api') and py.language_version().version_compare(
        '>=3.12',
    )
        limited_api_args = ['-DPy_LIMITED_API=0x030C0000']
    endif
    nanobind_lib = static_library(
        'nanobind',
        nb_staged[0],
        dependencies: py.dependency(),
        cpp_args: nb_compile + limited_api_args,
        pic: true,
    )
    nanobind_dep = declare_dependency(
        link_with: nanobind_lib,
        dependencies: py.dependency(),
        compile_args: nb_compile,
        link_args: nb_link,
    )
else
    nanobind_dep = dependency(
        'nanobind',
        default_options: [
            'python.allow_limited_api=' + get_option('python.allow_limited_api').to_string(),
        ],
    )
endif
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
        'src/pydseams/bond_network.py',
        'src/pydseams/features.py',
        'src/pydseams/frame.py',
        'src/pydseams/aseio.py',
        'src/pydseams/io.py',
        'src/pydseams/liquid_order.py',
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

# The sdist keeps wrap-git patches and the vesin and gpulite overlays.
meson.add_dist_script(py, files('scripts/dist_subprojects.py'))

# Shared-library paths. A build target has full_path(); an installed
# pkg-config dependency exposes the directory as libdir.
if host_machine.system() == 'darwin'
    _seams_soname = 'libseams.dylib'
    _lang_soname = 'libseams_lang.dylib'
elif host_machine.system() == 'windows'
    _seams_soname = 'seams.dll'
    _lang_soname = 'seams_lang.dll'
else
    _seams_soname = 'libseams.so'
    _lang_soname = 'libseams_lang.so'
endif
if use_installed_engine
    _seams_libdir = seams_dep.get_variable(pkgconfig: 'libdir')
    _seams_lib_path = _seams_libdir / _seams_soname
    _seams_lang_path = _seams_libdir / _lang_soname
    _built_libs = []
else
    _seams_libdir = meson.project_build_root() / 'subprojects' / 'seams-core' / 'src'
    _seams_lib_path = seams_api_lib.full_path()
    _seams_lang_path = seams_lang_lib.full_path()
    _built_libs = [seams_api_lib, seams_lang_lib]
endif

# These four use the shared libraries and synthetic cells. They do not
# read the seams-core checkout.
_lib_env = environment()
_lib_env.set('SEAMS_LIB', _seams_lib_path)
_lib_env.set('SEAMS_LANG_LIB', _seams_lang_path)
_lib_env.set('SEAMS_CORE_LIB', _seams_lib_path)
test(
    'python_chill_plus_ctypes',
    py,
    args: [files('tests/dseams/test_chill_plus_ctypes.py')],
    env: _lib_env,
    depends: _built_libs,
    timeout: 60,
)
test(
    'python_q4_q6_chill',
    py,
    args: [files('tests/dseams/test_q4_q6_chill.py')],
    env: _lib_env,
    depends: _built_libs,
    timeout: 120,
)
test(
    'python_frames_abi',
    py,
    args: [files('tests/dseams/test_frames_abi.py')],
    env: _lib_env,
    depends: _built_libs,
    timeout: 120,
)
test(
    'python_symmetry_ctypes',
    py,
    args: [files('tests/dseams/test_symmetry_ctypes.py')],
    env: _lib_env,
    depends: _built_libs,
    timeout: 120,
)

# Every module except the two that read the seams-core checkout.
_dseams_common_env = environment()
_dseams_common_env.set('DSEAMS_BUILD', meson.project_build_root())
_dseams_common_env.set('DSEAMS_SRC', meson.project_source_root())
_dseams_common_env.set('SEAMS_LIB', _seams_lib_path)
_dseams_common_env.set('SEAMS_LANG_LIB', _seams_lang_path)
test(
    'python_dseams_no_checkout',
    py,
    args: [
        files('tests/dseams/run_unittest.py'),
        '--exclude',
        'test_seams_package.py',
        'test_doc_snippets.py',
    ],
    env: _dseams_common_env,
    depends: [dseams_core] + _built_libs,
    timeout: 180,
)
test(
    'python_dseams_run_unittest_selftest',
    py,
    args: [files('tests/dseams/selftest_run_unittest.py')],
    env: _dseams_common_env,
    depends: [dseams_core] + _built_libs,
    timeout: 120,
)
# Twelve run_unittest.py processes at once for ten rounds, with
# MALLOC_PERTURB_ cycling through 0, 20, 139 and 255. They share
# <build>/python-stage, as the two tests above do when meson runs them side
# by side.
test(
    'python_dseams_stage_concurrency',
    py,
    args: [files('tests/dseams/stage_concurrency.py')],
    env: _dseams_common_env,
    depends: [dseams_core] + _built_libs,
    timeout: 600,
)

# Each of these needs the seams-core subproject checkout (input/, scripts/)
# or the libraries meson just built under that checkout.
if not use_installed_engine
    _engine_src = meson.project_source_root() / 'subprojects' / 'seams-core'
    _sub_libdir = meson.project_build_root() / 'subprojects' / 'seams-core' / 'src'

    # python_dseams: the modules that read the seams-core checkout.
    _dseams_env = environment()
    _dseams_env.set('DSEAMS_BUILD', meson.project_build_root())
    _dseams_env.set('DSEAMS_SRC', meson.project_source_root())
    _dseams_env.set('SEAMS_CORE_SRC', _engine_src)
    _dseams_env.set('SEAMS_LIB', _seams_lib_path)
    _dseams_env.set('SEAMS_LANG_LIB', _seams_lang_path)
    test(
        'python_dseams',
        py,
        args: [
            files('tests/dseams/run_unittest.py'),
            '--include',
            'test_seams_package.py',
            'test_doc_snippets.py',
        ],
        env: _dseams_env,
        depends: [dseams_core, seams_lang_lib],
        timeout: 180,
    )
    # python_seams_ctypes: reads input/traj and include/seams.h from the checkout.
    _ctypes_env = environment()
    _ctypes_env.set('PYTHONPATH', meson.project_source_root() / 'src')
    _ctypes_env.set('SEAMS_LIB', _seams_lib_path)
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

    # docs_api_python: ctypes reads input/ under the checkout.
    _examples_env = environment()
    _examples_env.set('SEAMS_CORE_SRC', _engine_src)
    _examples_env.set('SEAMS_LIBDIR', _sub_libdir)
    test(
        'docs_api_python',
        py,
        args: [files('examples/python/python_api.py'), _sub_libdir],
        env: _examples_env,
        depends: [seams_api_lib],
        timeout: 60,
    )

    # docs_minimal_python: genice_sI.lammpstrj and scripts/lammps_box.py.
    test(
        'docs_minimal_python',
        py,
        args: [
            files('examples/python/minimal_py.py'),
            _engine_src / 'input' / 'traj' / 'genice_sI.lammpstrj',
            _sub_libdir,
        ],
        env: _examples_env,
        depends: [seams_api_lib],
        timeout: 60,
    )

    # docs_core_counts: ice1h and genice_sII dumps plus scripts/ from the checkout.
    test(
        'docs_core_counts',
        py,
        args: [files('examples/python/core_counts.py')],
        env: _examples_env,
        depends: [seams_api_lib],
        timeout: 120,
    )

    # docs_lang_python: DSEAMS_ROOT/input/traj/ice1h.lammpstrj.
    _lang_tasks_env = environment()
    _lang_tasks_env.set('SEAMS_LANG_LIB', _seams_lang_path)
    _lang_tasks_env.set('SEAMS_CORE_LIB', _seams_lib_path)
    _lang_tasks_env.set('DSEAMS_ROOT', _engine_src)
    test(
        'docs_lang_python',
        py,
        args: [files('examples/python/lang_tasks_python.py')],
        env: _lang_tasks_env,
        depends: [seams_api_lib, seams_lang_lib],
        timeout: 120,
    )

    # docs_science_tutorials: the engine checkout and its build tree.
    # The C examples link against the engine with its compiler and OpenMP flags.
    _science_env = seams.get_variable('_snippet_env')
    if seams.get_variable('_chem_vendored')
        _loader_path = host_machine.system() == 'darwin' ? 'DYLD_LIBRARY_PATH' : 'LD_LIBRARY_PATH'
        _science_env.prepend(
            _loader_path,
            seams.get_variable('chemfiles_libdir'),
        )
    endif
    # The released engine's example driver accepts one compiler executable.
    _cc_cmd = meson.get_compiler('c').cmd_array()
    if _cc_cmd.length() > 1 and import('fs').name(_cc_cmd[0]) in [
        'sccache',
        'ccache',
    ]
        _science_env.set(
            seams.get_variable('_asan_wrap') ? 'REAL_CC' : 'CC',
            _cc_cmd[1],
        )
    endif
    _science_env.set('SEAMS_CORE_SRC', _engine_src)
    _science_env.set('SEAMS_LIBDIR', _sub_libdir)
    _science_env.set(
        'SEAMS_FIGURES',
        meson.project_build_root() / 'docs-science' / 'figures',
    )
    test(
        'docs_science_tutorials',
        py,
        args: [
            files('examples/python/science_tutorials.py'),
            _engine_src,
            meson.project_build_root() / 'subprojects' / 'seams-core',
        ],
        env: _science_env,
        depends: [seams_api_lib],
        timeout: 300,
    )
endif
