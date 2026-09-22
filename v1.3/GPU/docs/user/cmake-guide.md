# Tinker-HP CMake User Guide


## Supported GPU configurations

The native CMake path requires NVHPC through MPI wrappers and `nvcc`, uses the
GPU backend and the generic `2decomp_fft` backend.
It supports double precision and mixed precision. The public preset catalog has
34 configure presets: 32 static-core combinations and two constrained
shared-core presets. PLUMED, Colvars, and Deep-HP are separate feature choices
in this catalog; they are not combined in one public preset.

The common static Ninja workflow is:

```bash
cmake --preset gpu-double
cmake --build --preset gpu-double
```

The mixed-precision Fixed Precision Arithmetic (FPA) workflow is:

```bash
cmake --preset gpu-fpa
cmake --build --preset gpu-fpa
```

Feature variants use one suffix at a time. For example:

```bash
cmake --preset gpu-fpa-colvars
cmake --build --preset gpu-fpa-colvars

cmake --preset gpu-fpa-nn
cmake --build --preset gpu-fpa-nn

cmake --preset gpu-fpa-plumed
cmake --build --preset gpu-fpa-plumed
```

`gpu-fpa` produces executables such as `dynamic.fpa`; the Colvars and Deep-HP
variants add their corresponding executable markers.

For all-architecture compilation, use the explicit `all` token:

```bash
cmake --preset gpu-all-fpa
cmake --build --preset gpu-all-fpa
```

The same static choices are available with Unix Makefiles by adding `-make`:

```bash
cmake --preset gpu-fpa-make
cmake --build --preset gpu-fpa-make
```

Use the matching name for configuration and build. `cmake --list-presets`
shows the complete public catalog.

## Shared core linking

`TINKERHP_LINK_MODE` accepts `STATIC` (the default) and `SHARED`. Shared mode
builds only `tinkerhp_core` as a shared library; its implementation libraries
and the optional Colvars archive remain static and are compiled as PIC. The
installed executables use a relative runtime search path to find the installed
core library.

The public shared-core workflows are deliberately limited to the two
featureless, all-architecture Ninja presets:

```bash
cmake --preset gpu-all-double-shared
cmake --build --preset gpu-all-double-shared

cmake --preset gpu-all-fpa-shared
cmake --build --preset gpu-all-fpa-shared
```

The private SM75 presets select an explicit architecture for sandbox and
acceptance workflows. They are not part of the public 34-preset matrix. Run
the SM75 Colvars acceptance test after its build with:

```bash
ctest --preset gpu-sm75-fpa-shared-colvars
```

Shared-core preset availability is a configuration declaration, not evidence
that parallel MPI execution is safe. The current NVHPC compatibility evidence
defers shared-core MPI validation until NVHPC 26.1 or newer is available. Use
the maintained static SM75 presets for parallel development and validation.

## Preset naming at a glance

Public static names combine platform, architecture, precision, one optional
feature, link mode, and builder. Defaults are omitted: `gpu` is the platform,
native architecture, static link mode, and Ninja builder. `all`, `-make`, and
`-shared` are explicit when selected. The public precision token `fpa` means
mixed precision with FPA enabled. Feature suffixes are `-colvars`, `-nn`, and
`-plumed`; feature choices are mutually exclusive.

The complete matrix and mapping rules are maintained in the [CMake
architecture guide](../tech/cmake-arch.md). The matching Deep-HP workflow is:

```bash
cmake --preset gpu-fpa-nn
cmake --build --preset gpu-fpa-nn
```

It produces `dynamic_nn.fpa`. The public catalog does not define a combined
Deep-HP/Colvars preset; configure those features directly only where that
workflow is separately supported.

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
executables automatically. The selected interpreter's shared `libpython`
directory is recorded in the normal build and install runtime search paths for
both the plugin and NN executables; no plugin-specific `LD_LIBRARY_PATH`
setting is required. Installation deliberately does not copy Python
environments, packages, or models. The selected environment must remain at the
configured path (reconfigure and rebuild after moving or recreating it); an
ABI-compatible `LD_LIBRARY_PATH` override remains available for site
integration. Supply models through the existing `ML-MODEL` keyword. See the
[Deep-HP environment](../../deep-hp/README.md) for its lifecycle.

## Makefile Build
It remains possible to build using [Makefile](./makefile-build.md)
