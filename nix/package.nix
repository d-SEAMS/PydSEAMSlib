{
  lib,
  python3,
  meson,
  ninja,
  pkg-config,
  eigen,
  blas,
  lapack,
  libhwy,
  llvmPackages,
  stdenv,
  fetchFromGitHub,
  rustPlatform,
  seams-core-src,
}:

let
  geometryLibrary = { pname, version, owner, rev, hash }:
    rustPlatform.buildRustPackage rec {
      inherit pname version;
      src = fetchFromGitHub { inherit owner rev hash; repo = pname; };
      cargoLock.lockFile = "${src}/Cargo.lock";
      cargoBuildFlags = [ "--package" pname "--lib" ];
      cargoTestFlags = [ "--package" pname ];
      RAYON_NUM_THREADS = "8";
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
      '';
    };
  minimage = geometryLibrary {
    pname = "minimage";
    version = "0.1.4";
    owner = "lode-org";
    rev = "ce6dfb175f32459de1ee06dfeb70e8a718f61876";
    hash = "sha256-8Y9UFgj7hEiXcRNpUa0cHx6aT01QomJWCvvAuBfXya0=";
  };
  linkcell = geometryLibrary {
    pname = "linkcell";
    version = "0.3.11";
    owner = "d-SEAMS";
    rev = "2e26b83c0c9bbc07d09eb056011dd8e2cec6f4a5";
    hash = "sha256-36Zvbw/KgZw0VLUQstRIX27WXoan80JMKJNjkz9RgM0=";
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
      ../pyproject.toml
      ../README.md
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

  build-system = [ python3.pkgs.meson-python ];

  buildInputs = [
    eigen
    blas
    lapack
    libhwy
    minimage
    linkcell
    python3.pkgs.nanobind
  ]
  ++ lib.optionals stdenv.cc.isClang [ llvmPackages.openmp ];

  dependencies = [ python3.pkgs.numpy ];

  nativeCheckInputs = [
    python3.pkgs.pytest
    python3.pkgs.hypothesis
    python3.pkgs.ase
  ];

  # meson-python drives configure. The meson setup hook would leave
  # pypa looking at the build directory as if it were the project.
  dontUseMesonConfigure = true;
  dontUseMesonInstall = true;
  mesonAutoFeatures = "disabled";

  postPatch = ''
    mkdir -p subprojects
    cp -r ${seams-core-src} subprojects/seams-core
    chmod -R u+w subprojects/seams-core
    rm -rf subprojects/seams-core/subprojects

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
    description = "Python bindings for the d-SEAMS C++ engine";
    homepage = "https://github.com/d-SEAMS/PydSEAMSlib";
    license = lib.licenses.mit;
    platforms = lib.platforms.linux ++ lib.platforms.darwin;
  };
}
