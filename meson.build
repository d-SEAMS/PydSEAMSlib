project(
    'pydseams',
    'cpp',
    version: '2.11.1',
    default_options: [
        'warning_level=1',
        'cpp_std=c++20',
        'python.allow_limited_api=false',
    ],
    meson_version: '>=1.3.0',
)

# C++ engine. Bindings live here; seams-core is the library, not a monolith.
seams = subproject(
    'seams-core',
    default_options: [
        'default_library=static',
        'with_python=false',
        'with_linkcell=enabled',
        'with_minimage=enabled',
        'with_tests=false',
        'with_cli=false',
        'with_lua=disabled',
        'with_ira=auto',
    ],
)
yds_dep = seams.get_variable('yds_dep')

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
