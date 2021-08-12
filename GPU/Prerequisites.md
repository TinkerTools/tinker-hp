# Prerequisites

## Hardware
A relatively recent Nvidia GPU is mandatory for the GPU code. One with at least a compute capability above 5.0.

Nothing special is needed for the CPU code.


## Compilers, Operating Systems and Environment
   - A GNU or Intel FORTRAN compiler for the CPU code (Not recommended at this point)
   - The most recent [PGI compiler](https://www.pgroup.com/products/community.htm).
   - The nvidia CUDA C/C++ compiler for CUDA Files in GPU code
   - MPI Fortran wrapper on the selected compiler and compiled to be CUDA-Aware for device communications during multi-GPUs execution.
   - Linux or Windows 10 (Windows Subsystem for Linux) are preferred  
     __Note :__ PGI compiler does not support the GPU code on macOS.
   - [Nvidia HPC Package](https://developer.nvidia.com/nvidia-hpc-sdk-releases) Any version under __21.2 !__ 
     An issue involving "atomic operations" has been discovered with higher versions of the package (21.[3-7]). Luckily it is possible to bypass it by adding a compiler flag during the build step. Look for `add_options` variable configuration in `ci/install.sh` to be revert to the old behavior and fix the problem. However be aware that this is just a temporary fix. 
     This package contains everything you need to build Tinker-HP (GPU part). Previously listed items are already available within it. However we need to make sure the installation has been correctly done and the environment variables (PATH; CPATH & NVHPC) are correctly set and updated. Follow instructions described in `$PREFIX/modulefiles` with `$PREFIX` the installation directory if not.


## Mandatory Libraries
   - FFt libraries
      - `libfftw3` or `libfftw3f` in single or mixed precision mode.
        You will need to provide in that case FFTW install directory to both Tinker's and 2decompfft's Makefile.
        Exporting `FFTW=/path_to_fftw3` as an environment variable should be enough to make.
      - It is also possible to use the generic fft provided with your compiler (Default behavior).
   - Intel MKL Library (Only for host/CPU build!)
      - Export `MKLROOT=/path_to_mkl_library` inside your shell environment
      - `lmkl_intel_lp64 lmkl_sequential lmkl_core` are required
   - `2decomp_fft` library already provided in repository. Need to be compiled before Tinker
   - `libstdc++` Standard C++ library
      - Export `GNUROOT=/path_to_gnu`.
      This is not mandatory for a first and unique installation with `install.sh` but it is required for developers or when you need to recompile in an other configuration.
      e.g. `export GNUROOT=/usr` for conventional Linux systems
   - `CUDA` (Toolkit and libraries)
      - Cuda(C/C++) compiler `nvcc`
      - Recent CUDA libraries are required by the GPU code. Every CUDA version since 9.1 has been tested.
        They are delivered with PGI compiler.

### Notes
It is important to make sure that your environments (both build and run) are version consistent. In most cases, we have to pay attention our native CUDA version does not differ from PGI provided CUDA and/or NVIDIA driver installed. This matter will not happen with Nvidia HPC Package installed and loaded in your environment.
