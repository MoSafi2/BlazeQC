# BlazeQC

Fast quality control for sequencing data in Mojo.

## Installation

### From source (recommended)

Install [pixi](https://prefix.dev/docs/pixi/welcome) and [Mojo](https://docs.modular.com/mojo/), then:

```bash
git clone https://github.com/MoSafi2/BlazeQC.git
cd BlazeQC
pixi install
```

This creates an environment with Mojo and all dependencies (including BlazeSeq). Use `pixi run <task>` to run commands in that environment.

### From conda (after the package is published)

```bash
conda install -c conda-forge -c https://conda.modular.com/max blazeqc
```

### Mojo only (no pixi)

If you already have Mojo and dependencies installed:

```bash
mojo package blazeqc -o BlazeQC.mojopkg
```

Then use `-I BlazeQC.mojopkg` when building or running code that depends on BlazeQC.

## Build

### Build the Mojo package

```bash
pixi run build
```

This produces `BlazeQC.mojopkg` in the project root.

### Build the conda package

```bash
pixi run rbuild
```

This uses pixi-build and rattler-build to create a `.conda` package (e.g. under `.pixi/build/` or in the project root).

### Run tests

```bash
# Quick test (run test_helpers.mojo)
pixi run test

# All Mojo tests (build package first, then run test_*.mojo)
pixi run t

# Test the built conda package
pixi run rtest
```

## License

MIT
