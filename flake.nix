{
  description = "pydseams: Python bindings for the d-SEAMS C++ engine";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    # Source only: meson builds seams-core as a subproject. Same commit as
    # subprojects/seams-core.wrap; becomes github:d-SEAMS/seams-core/v2.14.0
    # once the tag exists.
    seams-core = {
      url = "github:HaoZeke/seams-core/628b0af94b4d3bf52426708a0d08436027c8e4fa";
      flake = false;
    };
  };

  outputs =
    { self, nixpkgs, seams-core }:
    let
      inherit (nixpkgs) lib;
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];
      forAllSystems = lib.genAttrs systems;
      pkgsFor = system: nixpkgs.legacyPackages.${system};
    in
    {
      packages = forAllSystems (
        system:
        let
          pkgs = pkgsFor system;
          pydseams = pkgs.callPackage ./nix/package.nix {
            python3 = pkgs.python312;
            seams-core-src = seams-core;
          };
        in
        {
          inherit pydseams;
          default = pydseams;
        }
      );

      checks = forAllSystems (system: {
        pydseams = self.packages.${system}.default;
        default = self.checks.${system}.pydseams;
      });

      devShells = forAllSystems (
        system:
        let
          pkgs = pkgsFor system;
        in
        {
          default = pkgs.mkShell {
            name = "pydseams-dev";
            inputsFrom = [ self.packages.${system}.default ];
            packages = with pkgs.python312.pkgs; [
              pytest
              hypothesis
            ];
          };
        }
      );

      formatter = forAllSystems (system: (pkgsFor system).nixfmt);
    };
}
