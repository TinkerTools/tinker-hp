# Prerequisites

## Hardware
A relatively recent Nvidia GPU is mandatory for the GPU code.
Nothing special is needed for the CPU code.


## Compilers, Operating Systems and Environnement
   - A GNU or Intel Fortran compiler for the CPU code (Not recommended at this point)
   - The most recent [PGI compiler](https://www.pgroup.com/products/community.htm).
   - The nvidia CUDA C/C++ compiler for CUDA Files in GPU code
   - MPI Fortran wrapper on the selected compiler
   - Linux or Windows 10 (Windows Subsystem for Linux) are preferred  
     PGI compiler does not support the GPU code on macOS.


## Libraries
   - FFT libraries
      - `libfftw3` or `libfftw3f` in single or mixed precision mode.
        You will need to provide in that case FFTW install directory to both Tinker's and 2decompfft's Makefile.
        Exporting `FFTW=/path_to_fftw3` as an environment variable should be sufficient to make.
      - It is also possible to use the generic fft provided with your compiler.
   - Intel MKL Library
      - Exporting `MKLROOT=/path_to_mkl_library` is sufficient enough
      - `lmkl_intel_lp64 lmkl_sequential lmkl_core` are required
   - `2decomp_fft` library already provided in repository
   - `libstdc++` Standard C++ library
      - Export `GNUROOT=/path_to_gnu`.
      This is not mandatory in first installation with `install.sh` but it is required for developers or when you need to recompile in an other configuration.  
      e.g. `export GNUROOT=/usr/lib64` for conventional linux systems
   - `CUDA` (Toolkit and libraries)
      - Cuda(C/C++) compiler `nvcc`
      - Recent CUDA libraries are required by the GPU code.
        They are delivered with PGI compiler.
