# Tinker-HP CMake Architecture

## Preset catalog

The public catalog contains exactly 34 configure presets:

- 32 static-core presets: two architecture choices (`native`, `all`), two
  precision choices (`double`, `fpa`), four mutually exclusive feature states
  (empty, `colvars`, `nn`, `plumed`), and two builders (Ninja, Unix Makefiles).
- Two shared-core exceptions: `gpu-all-double-shared` and
  `gpu-all-fpa-shared`. Both are featureless, all-architecture, Ninja presets.

The public name grammar is:

```text
platform-architecture-precision-features-link-builder
```

The `gpu` platform token is always present. Native architecture, static link
mode, the empty feature state, and Ninja are defaults and are omitted. The
`all` architecture is explicit, `shared` identifies shared core, and `-make`
identifies Unix Makefiles. The public precision token `fpa` maps to CMake
`TINKERHP_PRECISION=MIXED` with `TINKERHP_ENABLE_FPA=ON`; `double` maps to
`TINKERHP_PRECISION=DOUBLE`. Feature labels map to their CMake switches:

| Name token | Description label | CMake option |
| --- | --- | --- |
| `colvars` | `Colvars` | `TINKERHP_ENABLE_COLVARS=ON` |
| `nn` | `Deep-HP` | `TINKERHP_ENABLE_NN=ON` |
| `plumed` | `Plumed` | `TINKERHP_ENABLE_PLUMED=ON` |

Static link mode maps to `TINKERHP_LINK_MODE=STATIC`; shared mode maps to
`TINKERHP_LINK_MODE=SHARED`. Descriptions and `displayName` use the same form:

```text
GPU <Native Arch|All Arch> <double|fpa> [<feature label>] <Static Core|Shared Core> <Ninja|Makefile>
```

Every configure preset has a matching build preset and distinct binary
directory. Test presets remain representative rather than an exhaustive
matrix. Private SM75 presets remain in `CMakeUserPresets.json` and override
the public architecture with `75`; they do not expand the public catalog.

## Scope

The native build owns the maintained NVHPC/MPI-wrapper GPU path. The target
graph comprises static `tinkerhp_2decomp_fft` and `tinkerhp_wrapper`, Fortran
and CUDA object libraries, `tinkerhp_core`, and the executable targets. Source
inventory remains explicit in `cmake/TinkerHPSourceLists.cmake`, with
`source/Makefile.pgi` as the reference inventory.

## Link-mode contract

`TINKERHP_LINK_MODE=STATIC` retains the existing target graph. `SHARED` makes
only `tinkerhp_core` a shared library and compiles all objects and static
implementation libraries linked into it as PIC. Colvars remains its external
`libcolvars.a` Makefile target; shared mode adds `-fPIC` to that archive build.
The shared core is emitted and installed in the library destination, while
installed executables use a platform-relative sibling-library RPATH.

The public shared-core surface is intentionally narrow: only
`gpu-all-double-shared` and `gpu-all-fpa-shared` are exposed. Shared mode is
not opened to feature, Makefile, or native-architecture variants. Only
`tinkerhp_core` is shared; its implementation libraries and the optional
Colvars archive remain static PIC inputs.

This availability does not establish safe parallel MPI execution. NVHPC 24.7
has documented shared-object MPI named-COMMON corruption during two-rank
in-place collectives; shared-core MPI validation is deferred until NVHPC 26.1
or newer. See the [shared-core compatibility issue](../../.scratch/nvhpc-shared-core-mpi-compatibility/issues/01-upgrade-nvhpc-for-shared-core-mpi.md).

## Precision contract

`DOUBLE` uses `-r8` and the `.gpu` suffix. `MIXED` follows
`source/config.nvidia.mk`: it supplies `MIXED`, uses `-r4` for Fortran and CUDA
Fortran, appends `fastmath` to NVHPC GPU properties, and passes
`--use_fast_math` to CUDA C++. CUDA C++ also receives `USE_ERFC_HASTINGS`.

`TINKERHP_ENABLE_FPA` defaults to `ON`; in mixed precision it supplies
`USE_DETERMINISTIC_REDUCTION`. The default mixed suffix is `.fpa`; explicitly
disabling FPA changes it to `.mixed`. FPA is intentionally a no-op outside
mixed precision. `TINKERHP_ENABLE_ORTHO_BOX_ONLY` maps solely to
`ORTHOGONAL_BOX_SHAPE_ONLY`.

The suffix order places the Deep-HP marker before external-integration and
precision markers, so the public mixed/FPA Deep-HP preset produces `_nn.fpa`.
The public catalog does not define combined feature presets.

## External integrations

PLUMED is configured during CMake configuration and its external target remains
a dependency of `tinkerhp_core`. Colvars remains an external Makefile bridge.
Neither integration is refactored by the mixed/FPA option surface.

## Validation

CTest includes `tinkerhp-cmake-configuration-contract`, which configures the
public double, mixed/FPA, and mixed/non-FPA option combinations in the internal
configuration-only mode. It checks effective suffixes, shared-core PIC/install
behavior, Colvars archive handling, preset registration, and emitted compilation
definitions without claiming GPU numerical validation. The dedicated
`gpu-sm75-fpa-shared-colvars` CTest preset runs the built
`dynamic_clv.fpa` executable through one MPI rank from `examples/`.

## Deep-HP interface

`TINKERHP_ENABLE_NN` adds `NN_SUPPORT` to the existing Tinker-HP translation
units and configures a build-phase CFFI plugin target. Its inputs are copied to
the build tree, where `mlbuilder.py` emits a platform-native `mlinterface`
library. NN executables link it and receive build/install rpaths to the library
destination and to the selected interpreter's shared `libpython` directory; the
generated plugin receives the same Python runtime path. Installation includes
the plugin but excludes Python environments, packages, and models. The path
uses normal platform RUNPATH semantics, so an ABI-compatible loader override
is possible; moving the selected environment requires reconfiguration and
rebuild. Configuration fails when a shared `libpython` path cannot be
resolved.

Configuration resolves `deep-hp/.venv/bin/python` before an explicit
`DEEPHP_PYTHON_INTERPRETER`, creating the managed environment with locked `uv`
only when neither is selected. Managed Python must be 3.11; explicit Python may
be 3.10 or newer. The selected interpreter is checked only for required module
imports, so configuration does not initialize JAX devices or require a GPU.
The build report records its origin, interpreter, version, and managed CUDA
profile. This supersedes only the NN deferral in
`doc/spec/cmake-mixed-fpa-and-deep-hp.md`.
