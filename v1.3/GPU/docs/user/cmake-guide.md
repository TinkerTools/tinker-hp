# Tinker-HP CMake User Guide


## Supported GPU configurations

The native CMake path requires NVHPC through MPI wrappers and `nvcc`, uses the
GPU backend and the generic `2decomp_fft` backend.
It supports double precision and mixed precision. PLUMED and Colvars remain
external packages configured by CMake and built before Tinker-HP links.

Configure and build the baseline:

```bash
cmake --preset gpu-double
cmake --build --preset gpu-double
```

The mixed/FPA workflow is available as a preset:

```bash
cmake --preset gpu-fpa
cmake --build --preset gpu-fpa
```

The matching Colvars configuration is:

```bash
cmake --preset gpu-fpa-colvars
cmake --build --preset gpu-fpa-colvars
```

`gpu-fpa` produces executables such as `dynamic.fpa`; `gpu-fpa-colvars`
produces `dynamic_clv.fpa`.

## Shared core linking

`TINKERHP_LINK_MODE` accepts `STATIC` (the default) and `SHARED`. Shared mode
builds only `tinkerhp_core` as a shared library; its implementation libraries
and the optional Colvars archive remain static and are compiled as PIC. The
installed executables use a relative runtime search path to find the installed
core library.

The FPA shared-core workflows are:

```bash
cmake --preset gpu-fpa-shared
cmake --build --preset gpu-fpa-shared
cmake --preset gpu-fpa-shared-colvars
cmake --build --preset gpu-fpa-shared-colvars
```

`gpu-sm75-fpa-shared` and `gpu-sm75-fpa-shared-colvars` select SM75 explicitly
for environments where native GPU architecture detection is unavailable. Run
the SM75 Colvars acceptance test after its build with:

```bash
ctest --preset gpu-sm75-fpa-shared-colvars
```

The matching Deep-HP workflows are:

```bash
cmake --preset gpu-fpa-nn
cmake --build --preset gpu-fpa-nn
cmake --preset gpu-fpa-nn-colvars
cmake --build --preset gpu-fpa-nn-colvars
```

They produce `dynamic_nn.fpa` and `dynamic_nn_clv.fpa`, respectively.

## Precision and FPA

Set `TINKERHP_PRECISION=MIXED` to use the legacy Makefile-equivalent mixed
configuration. It uses the mixed preprocessor definition, single-precision
Fortran real handling, and GPU/CUDA fast-math settings.

Fixed Precision Arithmetic (FPA) is enabled by default and is meaningful only
for mixed precision. It maps to `USE_DETERMINISTIC_REDUCTION`. To compare with
non-FPA mixed precision, explicitly disable it:

```bash
cmake -S . -B build/gpu-mixed -G Ninja \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpic++ \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_CUDA_COMPILER=nvcc \
  -DTINKERHP_PRECISION=MIXED \
  -DTINKERHP_ENABLE_FPA=OFF
```

This configuration emits `.mixed`. FPA has no operational effect in double or
single precision.

## Orthogonal-box-only mode

`TINKERHP_ENABLE_ORTHO_BOX_ONLY=ON` adds only the legacy
`ORTHOGONAL_BOX_SHAPE_ONLY` definition. It does not alter compiler flags,
source selection, or link libraries.
Only Orthorombic boxes are available with this mode

## External packages

Colvars and PLUMED retain their existing external-target lifecycle. Enabling
PLUMED runs its `configure` step while CMake configures Tinker-HP; the
`tinkerhp_plumed` target builds and installs it before the Tinker-HP core links.
Colvars is similarly built by its existing Makefile bridge.

For PLUMED executables, expose its installed runtime library directory:

```bash
export LD_LIBRARY_PATH=<PLUMED_HOME>/lib:$LD_LIBRARY_PATH
```

## Build report

Each build tree contains `tinkerhp-build-report.txt`, including selected
precision, link mode, orthogonal-box mode, FPA setting, executable suffix,
compiler flags, and external integration details.

## Deep-HP environment

`TINKERHP_ENABLE_NN=ON` enables the maintained Deep-HP CFFI interface. CMake
selects an existing `deep-hp/.venv/bin/python` first, then an explicitly passed
`-DDEEPHP_PYTHON_INTERPRETER=/path/to/python`, and otherwise creates the local
environment with `uv sync --locked`. Managed environments require Python 3.11;
an explicit development interpreter requires Python 3.10 or newer. The selected
interpreter must import CFFI, NumPy, PyYAML, CuPy, JAX, and FeNNol.

Managed creation chooses its `cuda12` or `cuda13` profile from the CUDA runtime
reported by the NVHPC preflight. The generated native `mlinterface` plugin is
kept in the build tree, installed with Tinker-HP libraries, and found by NN
executables automatically. Installation deliberately does not copy Python
environments, packages, or models. Activate or otherwise expose the selected
environment before starting a simulation; supply models through the existing
`ML-MODEL` keyword. See the [Deep-HP environment](../../deep-hp/README.md) for
its lifecycle.

## Makefile Build
It remains possible to build using [Makefile](./makefile-build.md)
