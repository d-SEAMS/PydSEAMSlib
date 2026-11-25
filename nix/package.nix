{
  lib,
  python3,
  meson,
  ninja,
  pkg-config,
  cmake,
  eigen,
  blas,
  lapack,
  libhwy,
  llvmPackages,
  stdenv,
  fetchFromGitHub,
  fetchPypi,
  rustPlatform,
  seams-core-src,
}:

let
  # rustc passes -install_name only with -C rpath. Otherwise ld64 uses the
  # cargo -o path, target/<triple>/release/deps/lib….dylib, as LC_ID_DYLIB.
  # The next link copies that into LC_LOAD_DYLIB, and dyld opens it after
  # the build directory is gone. The Nix ld wrapper adds an rpath for
  # lib*.so only, so an @rpath id would not resolve these .dylib files.
  # Upstream Meson passes -Wl,-install_name for the cdylib; do the same,
  # aimed at the installed store path.
  #
  # Cargo invokes $RUSTC_WRAPPER $RUSTC_WORKSPACE_WRAPPER $RUSTC args.
  # `cargo auditable` is the workspace wrapper and treats its argv[1] as
  # the rustc path (the word rustc, or an absolute rustc). A -C flag in
  # that slot is executed as the compiler. The target-info probe also
  # lists --crate-type cdylib, so the flag is added only for a real
  # cdylib link and only after the rustc path.
  darwinCdylibHook = dylibFile:
    lib.optionalString stdenv.hostPlatform.isDarwin ''
      wrapdir=$(mktemp -d)
      cat > "$wrapdir/rustc-wrapper" <<'EOF'
      #!/bin/sh
      response_has() {
        kind=$1
        file=$2
        prev=
        while IFS= read -r line || [ -n "$line" ]; do
          if [ "$kind" = cdylib ]; then
            if [ "$prev" = "--crate-type" ] && [ "$line" = "cdylib" ]; then
              return 0
            fi
            case "$line" in
              --crate-type=cdylib) return 0 ;;
            esac
          else
            case "$line" in
              --print|--print=*) return 0 ;;
            esac
          fi
          prev=$line
        done < "$file"
        return 1
      }
      is_rustc() {
        [ "$1" = "rustc" ] || [ "$(basename "$1")" = "rustc" ]
      }
      next=$1
      shift
      rustc_tok=
      if ! is_rustc "$next"; then
        rustc_tok=$1
        shift
      fi
      cdylib=0
      probe=0
      prev=
      for arg in "$@"; do
        if [ "$prev" = "--crate-type" ] && [ "$arg" = "cdylib" ]; then
          cdylib=1
        fi
        case "$arg" in
          --crate-type=cdylib) cdylib=1 ;;
          --print|--print=*) probe=1 ;;
          @*)
            if response_has cdylib "''${arg#@}"; then
              cdylib=1
            fi
            if response_has print "''${arg#@}"; then
              probe=1
            fi
            ;;
        esac
        prev=$arg
      done
      if [ "$cdylib" = 1 ] && [ "$probe" = 0 ]; then
        set -- -C "link-arg=-Wl,-install_name,$DYLIB_INSTALL_NAME" "$@"
      fi
      if [ -n "$rustc_tok" ]; then
        exec "$next" "$rustc_tok" "$@"
      fi
      exec "$next" "$@"
      EOF
      chmod +x "$wrapdir/rustc-wrapper"
      export DYLIB_INSTALL_NAME="$out/lib/${dylibFile}"
      export RUSTC_WRAPPER="$wrapdir/rustc-wrapper"
    '';
  assertDarwinInstallName = dylibFile:
    lib.optionalString stdenv.hostPlatform.isDarwin ''
      if command -v otool >/dev/null 2>&1; then
        dump_id() { otool -D "$1"; }
      elif command -v llvm-otool >/dev/null 2>&1; then
        dump_id() { llvm-otool -D "$1"; }
      else
        echo "otool is required to check ${dylibFile}" >&2
        exit 1
      fi
      got=$(dump_id "$out/lib/${dylibFile}" | awk 'NR==2 { print $1 }')
      want="$out/lib/${dylibFile}"
      if [ "$got" != "$want" ]; then
        printf '%s: LC_ID_DYLIB is "%s", want "%s"\n' "${dylibFile}" "$got" "$want" >&2
        exit 1
      fi
    '';
  geometryLibrary = { pname, version, owner, rev, hash, cargoHash }:
    rustPlatform.buildRustPackage rec {
      inherit pname version cargoHash;
      src = fetchFromGitHub { inherit owner rev hash; repo = pname; };
      cargoBuildFlags = [ "--package" pname "--lib" ];
      cargoTestFlags = [ "--package" pname ];
      RAYON_NUM_THREADS = "8";
      LINKCELL_OPENMP = "0";
      preBuild = darwinCdylibHook "lib${pname}.dylib";
      postInstall = ''
        mkdir -p "$out/include" "$out/lib/pkgconfig"
        cp include/*.h include/*.hpp "$out/include/"
        cat > "$out/lib/pkgconfig/${pname}.pc" <<EOF
        prefix=$out
        libdir=$out/lib
        includedir=$out/include
        Name: ${pname}
        Description: Periodic geometry library
        Version: ${version}
        Libs: -L$out/lib -l${pname}
        Libs.private: -lpthread -lm ${lib.optionalString (!stdenv.hostPlatform.isDarwin) "-ldl"}
        Cflags: -I$out/include
        EOF
      '' + assertDarwinInstallName "lib${pname}.dylib";
    };
  minimage = geometryLibrary {
    pname = "minimage";
    version = "0.1.4";
    owner = "lode-org";
    rev = "ce6dfb175f32459de1ee06dfeb70e8a718f61876";
    hash = "sha256-8Y9UFgj7hEiXcRNpUa0cHx6aT01QomJWCvvAuBfXya0=";
    cargoHash = "sha256-8NSCv17ACkQbkBnycIZ3K9mUUOWwYX4u7OMii2F2znk=";
  };
  linkcell = geometryLibrary {
    pname = "linkcell";
    version = "0.6.0";
    owner = "d-SEAMS";
    rev = "v0.6.0";
    hash = "sha256-C1RPSOQHsBjfYxsw5Si6Q4+B+YGlQjNHitinndAv8XA=";
    cargoHash = "sha256-2QscNVUPmCrwO1cAHWnDZAtU8Hsr6Bl/WTO65WwDsC8=";
  };
  # seams-core reads .con frames through readcon-core. The binding
  # option with_readcon defaults to auto and finds this library by
  # pkg-config; the sandbox cannot run the cargo wrap.
  readcon-core = rustPlatform.buildRustPackage rec {
    pname = "readcon-core";
    version = "0.16.2";
    src = fetchFromGitHub {
      owner = "lode-org";
      repo = pname;
      rev = "9e12ea9727d5f3c85d8eeb4d927beb634e1bdb59";
      hash = "sha256-N3dH/eGnulSuq2aV0ghRul0v/OZZE5UBUUc/sUgnUB4=";
    };
    cargoHash = "sha256-VB5P5YC6JpXTy4cxAZc6kdiEUIw5WquMOlyfhn2NS3I=";
    cargoBuildFlags = [ "--lib" ];
    cargoTestFlags = [ "--lib" ];
    preBuild = darwinCdylibHook "libreadcon_core.dylib";
    postInstall = ''
      mkdir -p "$out/include" "$out/lib/pkgconfig"
      cp include/*.h include/*.hpp "$out/include/"
      cat > "$out/lib/pkgconfig/${pname}.pc" <<EOF
      prefix=$out
      libdir=$out/lib
      includedir=$out/include
      Name: ${pname}
      Description: CON and convel reader
      Version: ${version}
      Libs: -L$out/lib -lreadcon_core
      Libs.private: -lpthread -lm ${lib.optionalString (!stdenv.hostPlatform.isDarwin) "-ldl"}
      Cflags: -I$out/include
      EOF
    '' + assertDarwinInstallName "libreadcon_core.dylib";
  };
  minimagePython = python3.pkgs.buildPythonPackage rec {
    pname = "minimage";
    version = "0.2.3";
    pyproject = true;
    src = fetchPypi {
      inherit pname version;
      hash = "sha256-n5COD6ACMagZdmMAGKNJz9ID3k9GV681XkrM5VuvJM8=";
    };
    cargoDeps = rustPlatform.fetchCargoVendor {
      inherit pname version src;
      hash = "sha256-VIQGbCeIuxZTAO/3GMImJ6Ot1VWoJzwMCefGkG4zkaM=";
    };
    build-system = with rustPlatform; [ cargoSetupHook maturinBuildHook ];
    pythonImportsCheck = [ "minimage" ];
  };
  readconPython = python3.pkgs.buildPythonPackage rec {
    pname = "readcon-chemfiles";
    version = "0.23.2";
    pyproject = true;
    src = fetchPypi {
      pname = "readcon_chemfiles";
      inherit version;
      hash = "sha256-OOstinvqcBgu9XvjWe8EYuHJki9UnXeMH5LLMM4vrDc=";
    };
    cargoDeps = rustPlatform.fetchCargoVendor {
      inherit pname version src;
      hash = "sha256-elG4RKHcPFoXgw+UOBypvi49BfACTXiYaq17LnQnAH0=";
    };
    build-system = with rustPlatform; [ cargoSetupHook maturinBuildHook ];
    nativeBuildInputs = [ cmake ];
    dontUseCmakeConfigure = true;
    # readcon's C++ bridge scans Cargo's registry for the headers. Nix
    # stores those same crates under the fixed-output vendor directory.
    postPatch = ''
      substituteInPlace build.rs --replace-fail \
        'let mut roots = Vec::new();' \
        'let mut roots = vec![PathBuf::from("${cargoDeps}")];'
    '';
    pythonImportsCheck = [ "readcon" ];
  };
  nanobindSrc = fetchFromGitHub {
    owner = "wjakob";
    repo = "nanobind";
    rev = "v2.14.0";
    hash = "sha256-aa829i7/R5TN++/VVZDGrPFLBrfZcdXC/cfvovRX8/8=";
  };
  robinMapSrc = fetchFromGitHub {
    owner = "Tessil";
    repo = "robin-map";
    rev = "v1.4.0";
    hash = "sha256-Hkgxiq2i0TuqMK/bI5OMOn3LkmSE40NimDjK1FBZpsA=";
  };
in
python3.pkgs.buildPythonPackage {
  pname = "pydseamslib";
  version = (builtins.fromTOML (builtins.readFile ../pyproject.toml)).project.version;
  pyproject = true;

  src = lib.fileset.toSource {
    root = ./..;
    fileset = lib.fileset.unions [
      ../meson.build
      ../meson.options
      ../pyproject.toml
      ../README.md
      ../LICENSE
      ../examples
      ../docs/source/api.md
      ../scripts/stage_nanobind.py
      ../scripts/dist_subprojects.py
      ../scripts/yoda_public_names.py
      ../src
      ../tests
      ../subprojects/nanobind.wrap
      ../subprojects/robin-map.wrap
      ../subprojects/packagefiles
    ];
  };

  nativeBuildInputs = [
    meson
    ninja
    pkg-config
    python3.pkgs.meson-python
  ];

  build-system = [
    python3.pkgs.meson-python
    python3.pkgs.ninja
  ];

  buildInputs = [
    eigen
    blas
    lapack
    libhwy
    minimage
    linkcell
    readcon-core
    python3.pkgs.nanobind
  ]
  ++ lib.optionals stdenv.cc.isClang [ llvmPackages.openmp ];

  dependencies = [ python3.pkgs.numpy ];

  nativeCheckInputs = [
    python3.pkgs.pytestCheckHook
    python3.pkgs.hypothesis
    python3.pkgs.ase
    minimagePython
    readconPython
  ];

  # meson-python drives configure. The meson setup hook would leave
  # pypa looking at the build directory as if it were the project.
  dontUseMesonConfigure = true;
  dontUseMesonInstall = true;
  mesonAutoFeatures = "disabled";
  # The wheel default uses the Eigen wrap. The sandbox supplies Eigen
  # in buildInputs and resolves the native dependencies without downloads.
  pypaBuildFlags = [
    "--config-setting=setup-args=--force-fallback-for="
    "--config-setting=setup-args=--wrap-mode=nodownload"
  ];

  postPatch = ''
    mkdir -p subprojects
    cp -r ${seams-core-src} subprojects/seams-core
    chmod -R u+w subprojects/seams-core

    cp -r ${nanobindSrc} subprojects/nanobind-2.14.0
    chmod -R u+w subprojects/nanobind-2.14.0
    cp subprojects/packagefiles/nanobind/meson.build subprojects/nanobind-2.14.0/meson.build
    sed -i "/version: run_command/,/).stdout().strip(),/c\\  version: '2.14.0'," \
      subprojects/nanobind-2.14.0/meson.build

    cp -r ${robinMapSrc} subprojects/robin-map-1.4.0
    chmod -R u+w subprojects/robin-map-1.4.0
    cat > subprojects/robin-map-1.4.0/meson.build <<'EOF'
    project('robin-map', 'cpp', version: '1.4.0')
    robin_map_dep = declare_dependency(
      include_directories: include_directories('include'))
    meson.override_dependency('robin-map', robin_map_dep)
    meson.override_dependency('tsl-robin-map', robin_map_dep)
    EOF
  '';

  # meson-python already installs the extension; skip wrapping
  # seams-core as a second prefix.
  mesonInstallFlags = [ "--skip-subprojects" ];

  pytestFlags = [ "tests/python" ];

  pythonImportsCheck = [
    "pydseams"
    "pydseams.yoda"
    "pydseamslib"
  ];

  meta = {
    description = "Python bindings for the d-SEAMS engine: CHILL+, rings, cages, order parameters, and LAMMPS, XYZ, and CON readers";
    homepage = "https://github.com/d-SEAMS/PydSEAMSlib";
    license = lib.licenses.mit;
    platforms = lib.platforms.linux ++ lib.platforms.darwin;
  };
}
