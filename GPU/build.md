# Build Tinker-HP (GPU)

### Easy Build
A relatively easy way to build GPU version of Tinker-HP is to use the installation bash script. After setting your environment according to the Prerequisites, you can proceed to installation by typing in your shell.
```
$> pwd
#> /home/user/.../tinker-hp/GPU
$> ci/install.sh
```

As long as you stand for one configuration at the time, you can familiarize with the __source/Makefile.pgi's__ targets to customize your own build along.
In the event of an installation scipt failure, we recommend to follow the regular build procedure given below.

### Configuration options
   All those options are disabled by default.
   - `NVTX_SUPPORT` enables NVTX markers during profiling if set to 1. This option is not useful in release construct.
   - `NVSHMEM_SUPPORT` enables build with nvshmem library if set to 1. This option require `NVSHMEM_HOME` variable to be set on nvshmem install directory. It has not been tested with recent Nvidia HPC package.
   - `NO_MUTATION` disables soft-core computation if set to 1
   - `FPA_SUPPORT` enables fixed point arithmetic in non-Bonded Forces and energy reduction if set to 1. It requires mixed precision to be enabled through the Makefile precision variable.

You can force options at compile by editing in `source/Makefile.pgi:l62`. A proper way will be introduced in next developments.

### Using Makefile
You can almost do any construct in the code with `source/Makefile` linked to `source/Makefile.gpi`. Here are some useful tips to configure your build's type

##### configuration variables

  - `arch=(host|[device]) .or. (cpu|[gpu])` allows you to select target architecture
  - `prec=(single|mixed|[double]) .or. (s|m|[d])` changes the precision build
  - `prefix=(path/|[../bin])` controls the installation directory
  - `prog_suffix=(any_string|[])` append suffix to install binary. This variable may allow you to have multiple binaries in the installation directory.  
  For instance `analyze` or `analyze.gpu`
  - `compute_capability=(list|[60,70])` selects the device's compute capability for the construct. this variable shall accept a list of compute capability separated by a comma.   
  e.g. `compute_capability=35,60,75`
  - `cuda_version=([10.1])` contains the cuda version to be used for construct.  
    _PGI 19.10_ compilers only supports _cuda 9.2 10.0 and 10.1_ for OpenACC. To be consistent with the construct, we strongly recommend to use the same version of CUDA compiler as your OpenACC compiler.
  - `opt=(debug|debug+|[release])` decides the optimization level of the compilation

##### Some targets

  - `create_build`  
    Used In association with BUILD_DIR will create a new build directory
    with links to source files and _source/Makefile.pgi_. When BUILD_DIR is not specified, the default directory is `build/`. This target become useful when you want to keep object's files from a different build.  
    N.B. Depending on the precision you might have to rebuild the 22decomp_fft library.  
    e.g.  `make create_build BUIL_DIR=../build-mix`

  - `2decomp_fft_rebuild` ;  
    This phony target allow you to rebuild 2decomp_fft library in double precision
    Since the compilation uses the same directory for 2decomp_fft, you might want to
    be sure that your library is being compiled in the correct precision

  - `2decomp_fft_rebuild_single` ;  
    Same as previous target but excpet from the precision's construct which happen here to in single. There is only two precision modes for 2decomp_fft Library. Tinker mixed precision building requires 2decomp_fft single precision library.

  - `thrust_lib_rebuild` ;  
    [Re]Build the wrapper on CUDA thrust library

  - `all` ;  
    compile link and install `bin/dynamic` `bin/analyze` and `bin/minimize`. If you just want one program to be build just enter `make <program>`.  
    e.g. `make analyze dynamic`

  - `[dynamic|analyze|minimize].mixed` ;  
    compiles link and install analyze|dynamic|minimize.mixed. Those specific targets calls behind `make` on their main target. Only used one at the time to avoid over building.  
    e.g. `make analyze.mixed -j4` is fine

  - `[dynamic|analyze|minimize].single` ;  
    compiles link and install analyze|dynamic|minimize.single.

  - `[dynamic|analyze|minimize].cpu` ;  
    compiles link and install analyze|dynamic|minimize.cpu. CPU binaries

  - `libtinker` builds `libtinker.a`


### Custom Build
Let suggest that we want CPU binaries in double precision combined with an out-of-source build by hand. After linking `source/Makefile.pgi` to `source/Makefile`, we should follow the next script for an out-of-source build :
```
$> pwd
#> /home/user/PME
$> mkdir -p bin
$> cd source
$> make arch=cpu 2decomp_fft_rebuild          # build lib2decompfft.a
$> make thrust_lib_rebuild                    # build libwapper.a
$> make create_build BUILD_DIR=../build_CPU   # create build directory
$> cd ../build_CPU
$> make arch=cpu all -j6                      # build Tinker's program for CPU
$> ls ../bin                                  # list binaries
#> analyze dynamic minimize bar
```
or the following one for an in-source build
```
$> cd source
$> make arch=cpu prog_suffix=.cpu all -j4
$> ls ../bin
#> analyze.cpu dynamic.cpu minimize.cpu bar.cpu
```

You can also create a configuration bash file that store your configuration and run it with your desire targets.  
```
$> cat << EOF >Tconf.sh
   make FPA_SUPPORT=1 prec=mixed arch=device compute_capability=60 prog_suffix=.gmix $@
   EOF
$> chmod 740 Tconf.sh
$> cd source
$> ln -s ../Tconf.sh T
$> ./T 2decomp_fft_rebuild_single  # Rebuild 2decomp_fft library in single precision
$> ./T thrust_lib_rebuild
$> ./T all -j4    # Build binaries (*.gfix) for 60 compute_capability device using fixed precision
```
_______________________________
