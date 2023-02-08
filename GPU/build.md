# Build Tinker-HP (GPU)


## Direct Build
For most users, we recommend this relatively accessible method which lies on a bash installation script named `/tinker-hp/GPU/ci/install/sh`. After setting your environment according to the [Prerequisites](Prerequisites.md), change directory to the GPU folder, edit the configuration section of the installation script to select your build and proceed with the commands below.
```bash
$> pwd
#> /home/$user/.../tinker-hp/GPU
$> ./ci/install.sh
```
This script will create and install two types of binaries for every applications. Double and mixed precision program will be available by the end of the installation, in the installdir.   
As long as you stand for one configuration at the time, The installation script will serve as an ideal build assistant. In the event of the installation script failure, we recommend to follow one the procedures given below.

### Configuration
The configuration section of `ci/install.sh` lists some variables which are passed to the Makefile during the compiling stage. Edit this part according to your need, before executing the script. We keep this section updated as the development progresses. 
  - `dbuild `  for out-of-source build. If interested in an other configuration without losing the previous object files, you may use a different value like `$tinkerdir/build1`
  - `target_arch`  control the target platform architecture. Decide whether the build should target CPUs or GPUs
  - `c_c`   Stand for compute capability. It is specific to the GPU device. Default option target every Nvidia GPU from Pascal architecture to Ampere.
  - `build_plumed `  (switch) enables/disables *PLUMED*  interface. PLUMED is automatically built
  - `build_colvars ` (switch) enables/disables *COLVARS* interface. COLVAR is automatically built
  - `NN `            (switch) enables/disables  Neural Network python interface.  
  - `FPA `            (switch) enables/disables  fixed precision arithmetic in mixed precision (enabled by default)
  Requires a CUDA python environment consistent with the installation one before building.


## Using Makefile
This part concerns every developers or anyone familiar enough with source code compilation. `GPU/source/Makefile.pgi` allows every type of build. By editing configurations variables you can change default compilers, adapts flags, rename libraries, and so on. Here is the description of the useful variables need to configure your own build's.

#### configuration variables
  - `NVTX_SUPPORT` enables NVTX markers during profiling if set to 1. This option is not useful in release construct.
  - `NVSHMEM_SUPPORT` enables build with nvshmem library if set to 1. This option requires `NVSHMEM_HOME` variable to be set on nvshmem install directory. It has not been tested with recent Nvidia HPC package.
  - `NO_MUTATION` disables soft-core computation if set to 1
  - `FPA_SUPPORT` enables fixed point arithmetic in non-Bonded Forces and energy reduction if set to 1. It requires mixed precision to be enabled through the Makefile precision variable.
  - `PLUMED_SUPPORT` builds and link both Tinker-HP and Plumed Package supplied next to it.
  - `NN_SUPPORT`  builds, interfaces and links Tinker-HP with machine learning code
  - `arch=(host|[device]) .or. (cpu|[gpu])` allows you to select target architecture
  - `prec=(single|mixed|[double]) .or. (s|m|[d])` changes the precision build
  - `prefix=(path/|[../bin])` controls the installation directory
  - `prog_suffix=(any_string|[])` append suffix to install binary. This variable may allow you to have multiple binaries in the installation directory.  
  For instance `analyze` or `analyze.gpu`
  - `compute_capability=(list|[60,70])` selects the device's compute capability for the construct. this variable shall accept a list of compute capability separated by a comma.   
  e.g. `compute_capability=35,60,75,80`
  - `cuda_version=([10.1])` contains the cuda version to be used for construct.  
    _PGI 19.10_ compilers only supports _cuda 9.2 10.0 and 10.1_ for OpenACC. To be consistent with the construct, we strongly recommend to use the same version of CUDA compiler as your OpenACC compiler.
  - `opt=(debug|debug+|[release])` decides the optimization level of the compilation

#### Some targets

  - `create_build` ;  
    Used In association with BUILD_DIR will create a new build directory
    with links to source files and _source/Makefile.pgi_. When BUILD_DIR is not specified, the default directory is `build/`. This target become useful when you want to keep object's files from a different build.  
    **N.B.** Depending on the precision you might have to rebuild the `2decomp_fft` library.  
    e.g.  `make create_build BUIL_DIR=../build-mix`

  - `2decomp_fft_rebuild` ;  
    This phony target allow you to rebuild 2decomp_fft library in double precision
    Since the compilation uses the same directory for 2decomp_fft, you might want to
    be sure that your library is being compiled into the correct precision

  - `2decomp_fft_rebuild_single` ;  
    Same as previous target but excpet from the precision's construct which happen here to in single. There is only two precision modes for 2decomp_fft Library. Tinker mixed precision building requires 2decomp_fft single precision library.

  - `thrust` ;  
    Build the wrapper on CUDA thrust library

  - `thrust_lib_rebuild` ;  
    [Re]Build the wrapper on CUDA thrust library

  - `plumed` ;  
    Build and Install Plumed Libraries required by Tinker-HP

  - `libtinker` builds `libtinker.a`

  - `all` *(default target)* ;  
    compile link and install `bin/analyze` `bin/bar` `bin/dynamic` and `bin/minimize`. If you just want one program to be build just enter `make <program>`.  
    e.g. `make analyze dynamic`

  - `[dynamic|analyze|minimize].mixed` ;  *DEPRECATED*  
    compiles link and install analyze|dynamic|minimize.mixed *(mixed precision)*. Those specific targets call behind `make` on their main target. Only used one at the time to avoid over building.  
    e.g. `make analyze.mixed -j4` is fine

  - `[dynamic|analyze|minimize].single` ;  *DEPRECATED*   
    compiles link and install analyze|dynamic|minimize.single *(full single precision executables)*.

  - `[dynamic|analyze|minimize].cpu` ;  *DEPRECATED*  
    compiles link and install analyze|dynamic|minimize.cpu. CPU binaries


### Custom Build
Let suggest that we want CPU binaries in double precision combined with an out-of-source build by hand. After linking `source/Makefile.pgi` to `source/Makefile`, we should follow the next script for an out-of-source build :
```bash
$> pwd
#> /home/$user/.../tinker-hp/GPU
$> mkdir -p bin
$> cd source
$> make arch=cpu 2decomp_fft_rebuild          # build lib2decompfft.a
$> make thrust_lib_rebuild                    # build libwapper.a
$> make create_build BUILD_DIR=../build_CPU   # create build directory
$> cd ../build_CPU
$> make arch=cpu all -j6                      # build Tinker's program for CPU
$> ls ../bin                                  # list binaries
#> analyze dynamic minimize bar pimd
```
or the following one for an in-source build
```bash
$> cd source
$> make arch=cpu prog_suffix=.cpu all -j4
$> ls ../bin
#> analyze.cpu dynamic.cpu minimize.cpu bar.cpu
```

You can also create a configuration bash file that store your configuration and run it with your desire targets. 
```bash
$> cat << EOF >Tconf.sh
#!/bin/bash
make FPA_SUPPORT=1 COLVARS_SUPPORT=1 prec=mixed arch=device compute_capability=60 cuda_version=11.3 prog_suffix=.gmix $@
EOF
$> chmod 740 Tconf.sh
$> cd source
$> ln -s ../Tconf.sh T
$> ./T 2decomp_fft_rebuild_single  # Rebuild 2decomp_fft library in single precision
$> ./T -j8                         # Build binaries (*.gfix) for 70 compute_capability device using fixed precision
```
