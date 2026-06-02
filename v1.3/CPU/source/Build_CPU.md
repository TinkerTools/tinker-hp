# Prerequisites and Compilation for the CPU version

## Calculation Libraries

Tinker-HP requires the <span class="smallcaps">mkl</span> library, a FFT
library (such as <span class="smallcaps">fftw</span>) and a slightly
modified <span class="smallcaps">2decomp_fft</span> library (shipped
with Tinker-HP) in order to run. The
<span class="smallcaps">2decomp_fft</span> library enables parallel 3D
<span class="smallcaps">fft</span> computations based on 2d-pencils data
distribution (see
[<span class="smallcaps">2decomp_fft</span> site](http://www.2decomp.org))
based on a sequential implementation of
<span class="smallcaps">FFTs</span> such as the one provided by the
<span class="smallcaps">fftw</span> library.

## Parallel library

Tinker-HP also requires a recent enough
<span class="smallcaps">mpi</span> library supporting
<span class="smallcaps">mpi</span> 3.x standards such as non blocking
collectives. The code has been extensively tested with recent
Intel<span class="smallcaps">mpi</span> versions (such as intel
<span class="smallcaps">mpi</span> 5.1) and better performances have
been observed with this family of
<span class="smallcaps">mpi</span> implementation compared to other ones
such as Open<span class="smallcaps">mpi</span> .

# Installation

As Tinker-HP is shipped in source form, you need to compile it. This was
not always an easy task in the previous releases. Tinker-HP now uses a
`configure` script built with autotools packages from
<span class="smallcaps">Gnu</span> to ease the compilation and
installation process. The first thing you should do is to type :

`autoconf ; automake `

in the main directory. This will generate the `configure` file. Apart
from the usual options available with all `configure` scripts, there are
specific options for Tinker-HP.

    Usage: ./configure [OPTION]... [VAR=VALUE]...

    Optional Features:
      --enable-debug                Enable debug mode (check array bounds, implicit
                                    none, etc...). Should not be active in normal
                                    operations [default is no]
      --enable-skylake              Enable AVX512 Optimization for Skylake Processors
                                    [default is no]
      --enable-knl                  Enable AVX512 Optimization for KNL (Xeon Phi)
                                    Processors [default is no]
      --enable-fft-generic          Enable generic FFT mode [default is yes]
      --enable-fft-mkl              Enable MKL   FFT mode [default is no]
      --enable-fft-fftw3            Enable fftw3 FFT mode [default is no]
      --enable-fft-fftw3_f03        Enable fftw3_f03 FFT mode [default is no]
      --enable-plumed               Enable plumed interface [default is no]
      --enable-colvars              Enable Colvars interface [default is no]
    Optional Packages:
      --with-blaslib=<BLAS LIB>     Specify BLAS library [mkl, lapack or 
                                    /absolute/path/to/BLAS_library]
      --with-fftlib=<FFT LIB>       Specify a library for FFT called by 2decomp [mkl or
                                    fftw3 or /absolute/path/to/FFTW_library]

The ultimate goal of this script is to let you type

        ./configure ; make ; make install ; cd example ; ./ubiquitin2.run

and have everything compiled, installed and running.

## List of Options

As for all the `configure` scripts, you can choose the directory in
which the binaries will be copied. So, `configure` has `--prefix=<DIR>`.

Tinker-HP has a special interest to know if it will run on AVX-512
capable processors. So, `configure` has options for that:

- `--enable-slylake`

- `--enable-knl`

Recall that Tinker-HP needs to make a
<span class="smallcaps">fft</span> decomposition with a modified version
of the <span class="smallcaps">2decomp_fft</span> library, which in turn
needs a working <span class="smallcaps">fftw</span> library. This is why
you can find `configure` options about
<span class="smallcaps">fft</span> interface and library:

- `--enable-fft-generic`

- `--enable-fft-mkl`

- `--enable-fft-fftw3`

- `--enable-fft-fftw3_f03`

- `--with-fftlib=``<FFT LIB>`

Tinker-HP also needs some functions that resides in a working
<span class="smallcaps">blas</span> library. So, there is an option for
that:

- `--with-blaslib=``<BLAS LIB>`

Tinker-HP is now able to be interfaced with
<span class="smallcaps">Plumed</span>. So, there is an option for that :

- `--enable-plumed`

Tinker-HP is now able to be interfaced with
<span class="smallcaps">Colvars</span>. So, there is an option for that
:

- `--enable-colvars`

Finally, as there might be some execution problems, or compilation
problems for the users who develop code, Tinker-HP has an
`--enable-debug` option.

`configure` tries to find its path to reach a valid
<span class="smallcaps">mpi</span> compiler and a valid Fortran compiler
by unsetting the environment variables `$FC` and `$F77`, and reading the
environment variable `$PATH`. It also tries to find valid
<span class="smallcaps">fft</span> and
<span class="smallcaps">blas</span> libraries by reading `$FFTW` and
`$MKLROOT` or `$LAPACK` respectively. Most of the time, these
environment variables are defined through the `module` framework. As a
try, do a

        module available 2>&1 | less

to see if you have the `module` framework installed on your machine, and
to know what `module` you can load. `configure` then figures out how to
build the correct Makefiles.

By default, `configure` chooses the
<span class="smallcaps">mkl</span> library from Intel as the
<span class="smallcaps">blas</span> and
<span class="smallcaps">fftw3</span> libraries, sets the
`--enable-fft-mkl` option, does not make any processor optimization, and
disables debugging. Thus, typing `./``configure` give the same result as
if you have typed
`./``configure` `--enable-fft-mkl` `--with-blaslib=``mkl` `--with-fftlib=``mkl`.

## Using `configure`

If you want to have different settings than those used by default,
you’ll have to give `configure` more information. Be aware that
`configure` cannot magically guess anything. So, the information you
give must be precise and complete.

### Install Directory

By default, this is where you have unzipped and untarred the
distribution. If you want another place, use `--prefix=<DIR>`. You can
choose any directory you want, providing that you have permission to
create this directory and/or write in it.

### Processor optimization

The machine on which you compile is not always the one on which
Tinker-HP will run. If you know that Tinker-HP is going to run on
AVX-512 capable processors, you are strongly encouraged to use one of :

- for KNL processors (also known as Xeon-Phi)

- for Skylake processors.

as this will dramatically improve the execution speed. Otherwise, the
optimization will be done using the capabilities of the compilation
machine, as determined by the compiler.

### FFT interface

You can choose the interface of <span class="smallcaps">fftw</span> you
want to use. This has an effect on the
<span class="smallcaps">2decomp_fft</span> library. So :

- gives the generic <span class="smallcaps">fft</span>, with no call to
  <span class="smallcaps">fftw</span> library

- gives the
  <span class="smallcaps">mkl</span> <span class="smallcaps">fft</span>.
  It also automatically selects the
  <span class="smallcaps">mkl</span> library as the
  <span class="smallcaps">fftw</span> library.

- gives the fftw3 interface, and is designed to work with an external
  <span class="smallcaps">fftw</span> library.

- gives the fftw3 Fortran2003 interface, and is designed to work with an
  external <span class="smallcaps">fftw</span> library

### FFT Library

You can choose the <span class="smallcaps">fftw</span> library you want
to use. It can come from the <span class="smallcaps">mkl</span> suite,
or some <span class="smallcaps">fftw3</span> package (either system
installed, or compiled by you). So :

- :
  <div class="minipage">

  selects the <span class="smallcaps">mkl</span> library, but needs the
  variable `$MKLROOT` to be set to the absolute path of the
  <span class="smallcaps">mkl</span> library
  </div>

- :
  <div class="minipage">
  selects the <span class="smallcaps">mkl</span> library by giving the
  absolute path of the <span class="smallcaps">mkl</span> library
  </div>

- :
  <div class="minipage">
  selects the <span class="smallcaps">fftw3</span> library, but needs
  the variable `$FFTW` to be set to the absolute path of the
  <span class="smallcaps">fftw3</span> library
  </div>

- :
  <div class="minipage">
  selects the <span class="smallcaps">fftw3</span> library by giving the
  absolute path of the <span class="smallcaps">fftw3</span> library
  </div>

Here are typical commands you can type. If `$MKLROOT` has been correctly
set :

        ./configure --enable-fft-mkl --with-fftlib=mkl

If you wish to give the absolute path of the library :

        ./configure --enable-fft-mkl   --with-fftlib=/path/to/mkl/library
        ./configure --enable-fft-fftw3 --with-fftlib=/path/to/fftw3/library

These last commands can also be written this way :

        MKLROOT=/path/to/mkl/library ./configure --enable-fft-mkl   --with-fftlib=mkl
        FFTW=/path/to/fftw3/library  ./configure --enable-fft-fftw3 --with-fftlib=fftw3

### BLAS library

You can choose the <span class="smallcaps">blas</span> library you want
to use. It can come from the <span class="smallcaps">mkl</span> suite or
some <span class="smallcaps">lapack</span> package (either system
installed, or compiled by you). So :

- :
  <div class="minipage">
  selects the <span class="smallcaps">mkl</span> library, but needs the
  variable `$MKLROOT` to be set to the absolute path of the
  <span class="smallcaps">mkl</span> library
  </div>

- :
  <div class="minipage">
  selects the <span class="smallcaps">mkl</span> library, by giving the
  absolute path of the <span class="smallcaps">mkl</span> library
  </div>

- :
  <div class="minipage">
  selects the <span class="smallcaps">lapack</span> library, but needs
  the variable `$LAPACK`to be set to the absolute path of the
  <span class="smallcaps">lapack</span> library
  </div>

-  :
  <div class="minipage">
  selects the <span class="smallcaps">lapack</span> library, by giving
  the absolute path of the
  <span class="smallcaps">lapack</span> library.
  </div>

Here are typical commands you can type. If `$MKLROOT` has been correctly
set :

        ./configure --enable-fft-mkl --with-blas=mkl

If you wish to give the absolute path of the library :

        ./configure --enable-fft-mkl --with-blas=/path/to/mkl/library
        ./configure --enable-fft-mkl --with-blas=/path/to/lapack/library

These last commands can also be written this way :

        MKLROOT=/path/to/mkl/library    ./configure --enable-fft-mkl --with-blas=mkl
        LAPACK=/path/to/lapack/library  ./configure --enable-fft-mkl --with-blas=lapack

### DEBUG mode

This mode is primarily intended for developers, but can also be useful
if you experience errors while running Tinker-HP. Adding
`--enable-debug` to the `configure` command turns on **boundary
checking**, forces **implicit none**, sets the optimization level to
**0** (the lowest value) and enables **backtracing** and **all
warnings**. The compilation produces all the binaries and gives them the
`.debug` extension, so that you know that these binaries are not
optimized.

## Output of configure

`configure` produces a final log to resume what will be done. It
displays using colors (if available) all the information you gave, and
everything it has been able to catch from the environment.

Here is the result of a successful run of the `configure` command :

    FFTW=/usr/local/fftw-3.3.7/Intel/2018/impi/
    ./configure --enable-debug --enable-fft-fftw3 --enable-plumed --with-fftlib=fftw3
    --with-blaslib=lapack

where we give the absolute path of the
<span class="smallcaps">fftw3</span> library in the `$FFTW` variable,
ask for the DEBUG mode, enable the fftw3 interface, ask for the
<span class="smallcaps">Plumed</span> interface, use the
<span class="smallcaps">fftw3</span> library and want the lapack library
for <span class="smallcaps">blas</span>, assuming that the
`$LAPACK` variable is already set.

    configure:
    configure: **********************************************************************
    configure: **
    configure: ** Running Mode         : PLUMED DEBUG (binaries'extension is '.debug')
    configure: ** MPI Fortran Wrapper  : mpiifort
    configure: ** Fortran Compiler     : ifort
    configure: ** Fortran flags        : -fpp -DPLUMED -O0 -g -u -warn all -check bounds  
    configure: **                      : -no-ipo -no-prec-div -inline -heap-arrays 
    configure: **                      : -traceback
    configure: ** 2decomp Library      : -L ../2decomp_fft/src/ -l2decomp_fft
    configure: ** PLUMED  Library      : -L ../plumed/Intel/lib/ -lplumed -lplumedKernel
    configure: ** PLUMED  Includes     : -I ../plumed/Intel/include
    configure: ** FFTW3 Interface      : fftw3 of the FFTW3 library
    configure: ** FFTW3 Path           : /usr/local/fftw-3.3.7/Intel/2018/impi//lib
    configure: ** FFTW3 Includes       : -I /usr/local/fftw-3.3.7/Intel/2018/impi//include
    configure: ** FFTW3 Library        : -lfftw3
    configure: ** BLAS Type            : LAPACK
    configure: ** BLAS Path            : /usr/local/Libraries/lapack-3.8.0/Intel/2018
    configure: ** Prefix installation  : /home/lhj/neutron/Tinker/REL/PME/v1.2
    configure: ** Binaries location    : /home/lhj/neutron/Tinker/REL/PME/v1.2/bin
    configure: **
    configure: **********************************************************************
    configure:

This log confirms that we are in DEBUG mode and that we use the Intel
compiler `ifort` and the `mpiifort` wrapper from
Intel<span class="smallcaps">mpi</span>. As The
<span class="smallcaps">Plumed</span> interface has been selected, the
configure script shows the
<span class="smallcaps">Plumed</span> settings. The installation
directory where all binaries (with `.debug` extension) will be installed
is shown as well.

In the case of the <span class="smallcaps">Plumed</span> or the
<span class="smallcaps">Colvars</span> interface, the
`configure` command will configure the plumed or the Colvars library as
well.

## Making binaries

Once you are happy with the option you selected, it’s time to run the
`make` command, or even the `make install` command, which will compile
and install all at once.

As the compilation process takes care of the dependencies between
subroutines and modules, you can safely use the `-j` flag of the `make`
command to do parallel compilation. This would dramatically speedup the
compilation process.

Anyway, on modern machines, the compilation is not very long, except for
2 or 3 subroutines that can take up to 5(!) minutes to compile,
depending on the compiler you use and even on the fastest machines.
Everything should compile and link gracefully.

Using `make install` will copy the binaries into the directory you
selected with the `--prefix=<DIR>` option. The binaries produced will
have the following extension :

<div class="center">

| Mode             | Extension       |
|:-----------------|:----------------|
| `NORMAL`         | .prod           |
| `PLUMED NORMAL`  | \_plumed.prod   |
| `COLVARS NORMAL` | \_colvars.prod  |
| `DEBUG`          | .debug          |
| `PLUMED DEBUG`   | \_plumed.debug  |
| `COLVARS DEBUG`  | \_colvars.debug |

</div>

the `make install` command will also create 5 shell scripts, named after
the binaries, to setup the proper running environment. Don’t forget to
install the binaries you created, or you will not be able to run the
examples.

# Note for developers

## Writing new sources

We don’t want to impose you a unique style of writing. Indeed, we don’t
have one. But we just want to give you some rules we believe are
important for the consistency of Tinker-HP’s code.

#### File format

We use <span class="smallcaps">Fixed Form</span> format throughout all
the code, even though the code is written in
<span class="smallcaps">Fortran90</span>. This is mandatory. The
compilation process would not work otherwise.

#### Editing

We always use lowercase letters for code (except for printing purposes).
We indent all lines embedded in `do.....enddo`, `do.....while`, etc...
statements, or in `if...else...endif` constructs.

#### Variable declarations

You are required to use `implicit none`. If you compile in debug mode,
that will be enforced by the compiler.

The order we use to declare variables is:

1.  `integer` (4 bytes sized)

2.  `real` (8 bytes sized)

3.  `logical`

4.  `array` (in the same order)

5.  `character` (single string or array)

We always try not to mix different types of variables in the same
declaration line. This is not just because it is easier to read. That is
also because it is more memory efficient, particularly for
vectorization, where alignment in memory is crucial. `character`
variables should be put at the very end, since they can have arbitrary
lengths and almost never align to a memory boundary. We also try to
choose significant names for the variables.

#### Comments

We always begin a comment line by the `c` character. The `!` character
should only appear in the middle of a line. This is because the `!`
character at the beginning is reserved to introduce compiler directives.

If ever you create modules, please comment all the new variables you
create, like :

       c     maxvalue        atoms directly bonded to an atom
       c     maxgrp          user-defined groups of atoms
       c     maxtyp          force field atom type definitions
       c     maxclass        force field atom class definitions

In subroutines or functions, give as many comments as you believe is
needed to understand what your code is doing. That would be precious for
you, and for us as well.

## Compiling new sources

If you make development on Tinker-HP, it is likely that you would need
to add subroutines and modules in the source directory.

All modules should be put in files named `MOD_xxxxx.f`, even though they
are written in <span class="smallcaps">Fortran90</span>. All functions
and routines should be put in files with names beginning by a lowercase
letter and with `.f` extension. Please, try to find significant names (
`epolar1tcg2shortreal.f` is far better than `ep1tc2shre.f`). You are
required to follow this scheme as much as possible.

To compile your new sources, you should add them in the `Makefile.am`
file of the `source` directory. We’ve put some comments in this file, to
help you know where to put things. Search for the string `Add` in the
file.

There are 3 different cases[^1]:

1.  You created a new main program (like `analyze` or `dynamic`). Add a
    line

    `bin_PROGRAMS += yourmain`

    (with no extension) below the line `bin_PROGRAMS += testgrad`. Then,
    add the lines

    `yourmain_SOURCES = yourmain.f`

    and

    `yourmain_DEPENDENCIES = libtinkermod.a libtinkercalc.a `

    after the similar lines concerning `testgrad`.

    In the `Makefile.am` file of the `scripts` directory, you should
    also add the line :

    `-(cd $``bindir`` ; $(LN_S) -f $(LINK_TO) yourmain )`

    in the `install-exec-hook:` section, and add `yourmain` at the end
    of the `uninstall-binSCRIPTS:` section.

2.  You created new module(s). Add lines

    `libtinkermod_a_SOURCES += MOD_xxxxx.f`

    just below[^2] `libtinkermod_a_SOURCES += MOD_virial.f`

3.  You created functions and subroutines. Add lines like

    `libtinkercalc_a_SOURCES += yourexplicitfilename.f`

    just below[^3] the line `libtinkercalc_a_SOURCES += version.f`

You should now go in the main directory, where `configure.ac` resides,
and type `autoconf` and `automake`. `autoconf` should not generate any
message. `automake` will probably do, mainly because of a different
version than the one we used to create the distribution. In this case,
just type `aclocal` before running `automake` again. These 2 (or 3)
commands will generate a new `configure` script that takes care of your
new sources. You should then run `./``configure`[^4], compile, install
and enjoy debugging your code.

# Executables

After having successfully compiled the code, five executable files
should be present in the install directory: `analyze`, `bar`, `dynamic`,
`testgrad` and `minimize`[^5], which are the analogous of the binaries
of the Tinker-8.4 release and require similarly a geometry (given by a
\*.xyz file), a simulation setup (given by a \*.key file) and possibly a
restart (given by a \*.dyn file) for the `dynamic` program.

All these executables must run in the same environment you had during
the compilation phase. That means the same set of modules, or the
correct `LIBRARY_PATH`. They should be launched with the
`mpirun -np `` x` prefix in order to run in parallel with `x`
<span class="smallcaps">mpi</span> processes.

## General remarks

The only boundary conditions that are available in this release are
periodic boundary conditions treated with Particle Mesh Ewald.
General triclinic unit cells can be used.

Classical force fields such as AMBER, CHARMM and OPLS are available in
Tinker-HP as well as polarizable force fields such as AMOEBA, AMOEBA+ and HIPPO.

## `analyze` 

The `analyze` executable allows potential energy analysis. Compared to
the Tinker-8.4 software, the only option compatible with this binary is
"e".

For example the command line:

- `mpirun -np`` 16 ../bin/analyze dhfr2 e`\
  will give you as an output the potential energy terms of the geometry
  given by a dhfr2.xyz file and with the simulation setup given by the
  dhfr2.key file. Furthermore, this computation will run on 16
  <span class="smallcaps">mpi</span> processes.

## `testgrad` 

The `testgrad` program is absolutely equivalent to the one of the
Tinker-8.4 release: it allows the output of the components of the
analytical and/or numerical gradients of the different energy terms.

For example, the command line:

- `mpirun -np`` 16 ../bin/testgrad dhfr2 Y Y 0.0001 Y`\
  will give you as an output all the analytical and numerical gradients
  (computed with an increment of 0.0001 Angstroms for the positions of
  the atoms) of all the energy terms of the dhfr2 system.

## `minimize`

The `minimize` program computes energy minimization starting from a
given structure, using a low memory quasi-newton BFGS algorithm as in
Tinker-8.4. The command line used should give the numerical threshold
for the convergence of the algorithm.

For example, the command line:

- `mpirun -np`` 16 ../bin/minimize dhfr2 0.1`\
  will compute energy minimization on the dhfr2 structure until the RMS
  on the gradient is inferior to 0.1. The new geometry will be written
  at each iteration of the algorithm in the file dhfr2.xyz_2.


# Examples

4 examples of systems with associated \*.key files are given in the
distribution: ubiquitin2, dhfr2, puddle and pond. The sizes of these
systems are respectively: 9737, 23558, 96000 and 288000 atoms, making
them good various benchmarks for the program.

4 different setups are given for the ubiquitin system with 4 different
key files:

|  |  |  |  |
|:---|:---|:---|:---|
| 1\) | **ubiquitin2.key** | : | regular (langevin with BAOAB integration based) 2 *fs* respa (bonded/non-bonded split) computations with DC-JI/DIIS as a polarization solver |
| 2\) | **ubiquitin2tcg.key** | : | 2 *fs* respa computations with TCG2 (with a diagonal preconditioner, no guess and a peek step with $`\omega=1`$ as a polarization solver |
| 3\) | **ubiquitin2respa1.key** | : | 6 *fs* respa1 (bonded/short range non-bonded/long range non-bonded split) langevin with BAOAB integration computations with DC-JI/DIIS as a short and total polarization solver |
| 4\) | **ubiquitin2respa1tcg.key** | : | 10 *fs* respa1 langevin with BAOAB integration computations with heavy hydrogen, TCG1 (with a diagonal preconditioner, no guess and no peek step) as a short range polarization solver and DC-JI/DIIS as a total polarization solver |

A fifth example is reserved for debug purposes. It’s exactly the same as
the first one. `./ubiquitin2.debug.run` runs the `dynamic`.`debug`
binary.

# Support

Tinker-HP is maintained by few people. That means we cannot promise you
to answer in a minute to your requests. Anyway, if you have any question
or need any support for Tinker-HP, feel free to send a mail to our team
at `TinkerHP_Support@ip2ct.upmc.fr`. We will answer as soon as we can,
providing that we can !

[^1]: Of course, you can match all three at the same time!

[^2]: As the compilation process takes care of all the dependencies, the
    positions of the lines you add are not really significant. But
    putting the new lines at the end is just a way of remembering they
    are – well – new.

[^3]: Same remark as above.

[^4]: Presumably with the `--enable-debug` option flag.

[^5]: If you have ever compiled with `--enable-debug` before, you should
    have 5 more binaries.
