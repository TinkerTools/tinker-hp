# Prerequisites


## Hardware
A relatively recent Nvidia GPU is mandatory for the GPU code. One with at least a compute capability above 5.0.

Nothing special is needed for the CPU code.


## Compilers, Operating Systems
   - A GNU or Intel FORTRAN compiler for the CPU code (Not recommended at this point)
   - The most recent [PGI compiler](https://www.pgroup.com/products/community.htm).
   - The nvidia CUDA C/C++ compiler for CUDA Files in GPU code
   - MPI Fortran wrapper on the selected compiler and compiled to be CUDA-Aware for device communications during multi-GPUs execution.
   - Linux or Windows 10 (Windows Subsystem for Linux) are preferred  
     __Note :__ PGI compiler does not support the GPU code on macOS.
   - [Nvidia HPC Package](https://developer.nvidia.com/nvidia-hpc-sdk-releases) Any version under __22.7 !__ is operational.  
     An issue involving "atomic operations" has been discovered with versions of the package (21.[3-9]). Luckily it is possible to bypass it by adding a compiler flag during the build step. Look for `add_options_f` variable configuration inside `ci/install.sh` to be revert to the old behavior and fix the problem. However be aware that this is just a temporary fix.  
     Starting from version 22.XX, the previously described issue has been solved by the developers, and the OpenACC compiler bug lately detected within 22.XX version has also been fixed. It involves an ambiguous OpenACC compilation when calling an acc device subroutine compiled in an other translation unit.  
     This package contains everything you need to build Tinker-HP (GPU code). Previously listed items are already available within it. However, we need to make sure the installation has been correctly done and the environment variables (PATH; CPATH & NVHPC) are set and updated. Follow instructions described in `$PREFIX/modulefiles` with `$PREFIX` the installation directory of your package.

## Issues
   #### SDK version
   Lately, we have been reporting some unexpected behaviors on binaries generated with Compiler versions strictly above 22.7. Depending on the machine, some programs may or may not crash at start. Discussions are initiated in order to identify and resolve the issue. Meanwhile, here are some stable setup able to run Tinker-HP

   - HPC-SDK 22.7 + cuda11.7 + GNU-11.x.x
   - HPC-SDK 22.7 + cuda11.0 + GNU-9.x.x
   - HPC-SDK 22.2 + cuda11.6 + GNU-9.x.x
   - HPC-SDK 22.2 + cuda11.0 + GNU-9.x.x
   - HPC-SDK 21.9 + cuda11.4 + GNU-8.x.x
   - HPC-SDK 21.2 + cuda11.2 + GNU-8.x.x

   If you have a multi-CUDA HPC-SDK installation, you can reconfigure the compilers to use the desired version of CUDA by editing/creating the `localrc` file inside the `nvfortran` compiler directory. Fill or complete the file with the instruction.
   ```
   set CUDAVERSION=11.x;
   ```
   In the case where you do not have write access in the directory, you can also edit a local file and _export_ an env variable named `NVLOCALRC` to hold the file path.


## Python Environment 
If you plan to combine machine learning with MD, we provide an environment file `tinker-hp/GPU/tinkerml.yaml` which installs python modules for machine learning. This specific python environment requires [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) to be installed on the host machine. Install the environment according to the next instructions before proceeding to Tinker-HP.

To create the environment, run in your shell: `conda env create -f tinkerml.yaml`; <br />
Activate or deactivate it respectively with `conda activate tinkerml` or `conda deactivate`;  <br />

Composition of the environment:  
* [Pytorch](https://pytorch.org/), [TensorFlow](https://www.tensorflow.org/) are the building block of most of machine learning potential libraries.
* Deepmd-kit and libdeepmd are used for [DeePMD](https://docs.deepmodeling.com/projects/deepmd/en/master/index.html) models.
* Our [TorchANI](https://aiqm.github.io/torchani/)-based library composed of lot of new features, more coming soon!


## Mandatory Libraries
   - FFt libraries
      - `libfftw3` or `libfftw3f` in single or mixed precision mode.
        You will need to provide in that case FFTW install directory to both Tinker's and 2decompfft's Makefile.
        Setting `FFTW=/path_to_fftw3` as an environment variable should be enough to make.
      - It is also possible to use the generic fft provided with your compiler (Default behavior).
   - Intel MKL Library (Only for host/CPU build!)
      - Set `MKLROOT=/path_to_mkl_library` inside your shell environment
      - `libmkl_intel_lp64 libmkl_sequential libmkl_core` are required
   - `2decomp_fft` library already provided inside the repository.
   - Standard C++ library
      - Set `GNUROOT=/path_to_gnu`. This is not mandatory for a first and unique installation with `ci/install.sh` but it is required for developers or when you need to recompile in an other configuration.  
      e.g. `export GNUROOT=/usr` for conventional Linux systems. The aim is to fetch `libstdc++` from `${GNUROOT}/lib64`
   - `CUDA` (Toolkit and libraries)
      - Cuda(C/C++) compiler `nvcc`
      - Recent CUDA libraries are required by the GPU code. Every CUDA version since 9.1 has been tested.
        From now on, CUDA Toolkit is a part of NVIDIA HPC-SDK package.

:warning: ** !! Warning !!**  
If you are building with machine learning interface, please ensure the consistency of your environment between the one of Tinker-HP (HPC-SDK) and python ; especially the installed version of CUDA. Default python environment comes with cuda11.3 version `(tinkerml.yaml)`. Therefore, it is mandatory to use an NVIDIA HPC-SDK package which contains cuda11.3?  
Should that instruction being discarded, unexpected behaviors may occur while running Tinker-HP.


## Notes
It is necessary to make sure that your environments (both build and run) are version consistent. In most cases, we have to pay attention, our native CUDA version does not differ from PGI provided CUDA and, is compatible with NVIDIA installed driver. Naturally, this issue will not occur with Nvidia HPC-SDK Package.
