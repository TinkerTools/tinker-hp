# Prerequisites

This page lists the packages and environment requirements for the maintained
GPU CMake path. Configuration and build instructions are in the
[CMake user guide](cmake-guide.md).

## Platform and hardware

Use a Linux system with an NVIDIA GPU and an NVIDIA driver compatible with the
CUDA runtime selected below. WSL2 is suitable when the installed NVIDIA driver,
NVIDIA HPC SDK, and MPI implementation support the selected WSL2 environment.
macOS and the CPU-only build are not supported by the maintained CMake path.

## Required packages

| Package | Requirement | Purpose |
| --- | --- | --- |
| NVIDIA HPC SDK (NVHPC) | `nvfortran` and `nvc++` | GPU Fortran/C++ compilation and OpenACC support. |
| CUDA Toolkit from the same toolchain | `nvcc` | CUDA C/C++ compilation and CUDA libraries. |
| MPI implementation | CUDA-aware; wrappers named `mpif90`, `mpic++`/`mpicxx`, and `mpicc`; compiled with NVHPC | Parallel execution and device-to-device MPI communication. CUDA awareness is required for multi-GPU runs. |
| CMake | version 3.31 or newer | Native project configuration. |
| Ninja | current version | CMake generator used by the supplied presets. |
| BLAS and LAPACK development libraries | system-provided implementation | Required by the CMake project. OpenBLAS plus LAPACK is a suitable choice. |
| GNU compiler collection | available while installing NVHPC | Required by NVIDIA's Linux HPC SDK installer. It is not the compiler used for the maintained Tinker-HP GPU path. |

The repository already contains `2decomp_fft`; do not install it separately.
The generic FFT backend is used by the maintained CMake path, so FFTW and MKL
are not prerequisites.

### Quick installation route

On a Debian/Ubuntu workstation, install the general build tools and numerical
libraries with:

```bash
sudo apt update
sudo apt install build-essential cmake ninja-build libopenblas-dev liblapack-dev
```

Check that the distribution's CMake is at least 3.31. If it is older, install a
current binary release from the [CMake download page](https://cmake.org/download/)
or use the CMake package provided by your site.

Install the NVIDIA HPC SDK using the packages or tar installer from the
[NVIDIA HPC SDK download page](https://developer.nvidia.com/hpc-sdk-downloads),
then load its environment module or source its environment script. The NVIDIA
[installation guide](https://docs.nvidia.com/hpc-sdk/installation-guide/index.html)
describes both local and network installations. On a cluster, the simplest
route is normally to load the site-provided NVHPC and NVHPC-built MPI modules.

Install or select an MPI distribution built against that NVHPC installation.
A distribution MPI package compiled with GCC is not a substitute: its wrappers
must invoke `nvfortran` and `nvc++`. Ask the system administrator for the
matching CUDA-aware MPI module if one is not already available.

After loading the toolchain, these checks should resolve to the same NVHPC and
MPI installation:

```bash
cmake --version
ninja --version
nvfortran --version
nvc++ --version
nvcc --version
mpif90 --version
mpic++ --version
```

## Deep-HP additions

Deep-HP is optional. When it is enabled, install [uv](https://docs.astral.sh/uv/getting-started/installation/)
in addition to the packages above. On Linux, its official installer is:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Ensure the installation directory is on `PATH` and confirm it with `uv --version`.
The managed Deep-HP environment requires Python 3.11; `uv` can install the
required Python automatically. Do not manually install the Python libraries:
the repository's locked `deep-hp` environment supplies CFFI, NumPy, PyYAML,
CuPy, JAX/JAXlib, and FeNNol. It selects either its CUDA 12 or CUDA 13
dependency profile from the NVHPC CUDA runtime.

An explicitly supplied Deep-HP development interpreter may use Python 3.10 or
newer, but it must already provide CFFI, NumPy, PyYAML, CuPy, JAX, and FeNNol.

## Compatibility requirements and known exceptions

- NVHPC 22.9 through 24.5 is blocked by the maintained CMake path because that
  range is known to be unstable for this project.
- NVHPC 21.3 through 21.9 needs the `-Mnollvm` workaround. CMake enables it by
  default; disabling the workaround is unsupported unless the issue is being
  investigated deliberately.
- The CUDA runtime selected by NVHPC's Fortran compiler and the `nvcc` compiler
  must have the exact same major.minor version. CMake stops during configuration
  when they differ; do not mix a system `nvcc` with a different CUDA version
  from the HPC SDK.
- For Deep-HP's managed environment, the detected NVHPC CUDA major version must
  be 12 or 13. Its selected `cuda12` or `cuda13` Python profile must match that
  major version.
- The installed NVIDIA driver must support the CUDA runtime chosen by NVHPC.
- For multi-GPU jobs, MPI must be CUDA-aware and must use the same NVHPC
  toolchain. A non-CUDA-aware or GCC-built MPI installation can compile or link
  incorrectly and cannot provide the required device communication path.

If one of these checks fails, correct the package selection or module order
before following the [CMake user guide](cmake-guide.md).
