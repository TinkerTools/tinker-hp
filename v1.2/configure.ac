#
#     Sorbonne University
#     Washington University in Saint Louis
#     University of Texas at Austin
#

#
#  configure.ac file to generate configure and Makefiles for Tinker-HP
#

AC_PREREQ([2.69])
AC_INIT([Tinker-HP], [1.2], [tinkersupp@ip2ct.upmc.fr])
#AM_INIT_AUTOMAKE(subdir-objects)
AM_INIT_AUTOMAKE()

#
# Unset FC and F77 environment variable
# Some modules set these variables, but that prevents configure from finding
# the MPI compiler

unset FC
unset F77

#
# Find the default ranlib program (help linking static libraries)
#

AC_PROG_CXX([mpiicpc mpic++])
AC_PROG_RANLIB

#
# Set the default prefix to the directory where configure runs
#
AC_PREFIX_DEFAULT($PWD)

# 
# Grab the prefix given by the --prefix option of configure
# If none, default to the... default !
# That's just for the NOTICE
#

PREFIX=$prefix
AS_IF([test "x$prefix" = xNONE],
      [PREFIX=$ac_default_prefix],
      [PREFIX=$prefix]
     )

#
# We are using Fortran 
#
AC_LANG(Fortran)

#
# No color
#

BOLD=""
BLACKBACK=""
RED=""
GRN=""
YEL=""
CLR=""
PADD="                                                                                                                       "
#
# Check for tput to get colors
#

AC_CHECK_PROG(tput,tput,"ok","no")

#
# Define some colors  if tput is available
#

AS_IF([test $tput = ok],
      [BOLD="$(tput bold)"
       BLACKBACK="$(tput setab 0)"
       RED="$BOLD$BLACKBACK$(tput setaf 1)"
       GRN="$BOLD$BLACKBACK$(tput setaf 2)"
       YEL="$BOLD$BLACKBACK$(tput setaf 3)"
       CLR=$(tput sgr0)
      ]
     )

#
# Define somme messages
#

FORT_COMP_FAILMSG="${RED}No known compiler found !!! Must have ${GRN}ifort${RED}, "
FORT_COMP_FAILMSG="$FORT_COMP_FAILMSG ${GRN}pgfortran ${RED} or ${GRN}gfortran$CLR"
CHECKBOUNDS_FAILMSG="${RED} No Flag available to check array bounds !!! $CLR"
NO_FFT_FAILMSG="${RED}Must have a FFT Library if FFTW_TYPE is not generic !$CLR"
PYTHON_FAILMSG="${RED}No python found !!! Must have it to compile sources ! $CLR"

#
# Check for python2 to run Tinker_deps.py
#

AC_CHECK_PROGS([PYTHON],[python  python2  python3],["no"])

AS_IF([test "x$PYTHON" = "xno"],
      [AC_MSG_FAILURE([$PYTHON_FAILMSG])],
      [:]
     )

#
# Find suitable MPI compiler for .f90 and .f sources
#

AC_PROG_LN_S                                       # Check if ln is available
AC_PROG_FC([mpiifort mpif90 ])
AC_PROG_F77([mpiifort mpif90 ])
AC_MSG_NOTICE(["$YEL$FC$GRN wrapper will be used to compile Tinker-HP${CLR}"])
MPIFC=$FC

#
# Checks for the Real Fortran Compiler (mpiifort and mpif90 are wrappers)
#

AC_CHECK_PROGS([REALFC],[ifort pgfortran gfortran],["none"])
AS_IF([test "x$REALFC" = "xnone" ],
      [AC_MSG_FAILURE([$FORT_COMP_FAILMSG])],
      [AC_MSG_NOTICE(["${GRN}The wrapped compiler is $YEL$REALFC$CLR"])]
     )
 
#
# Find the flag to set DEFINES for the preprocessing of fortran sources

AC_FC_PP_DEFINE

#
# All sources will be in fixed forms
#

AC_FC_FIXEDFORM

#
# Define the 2decomp_fft location
# and Include and Load flags

LIB2DECOMP_PATH="../2decomp_fft/src/"
LIB2DECOMP_INC="-I $LIB2DECOMP_PATH"
LIB2DECOMP_LIB="-L $LIB2DECOMP_PATH -l2decomp_fft"

#
#Define some variables 
#

#
# compilation flags for intel, portland and gnu compilers
#

NVMPICXX=""
FORT_FLAGSIF="-no-ipo -no-prec-div -inline -heap-arrays -traceback" 
FORT_FLAGSPF="-traceback -Mbackslash -Mdalign -gopt" 
FORT_FLAGSGF="-C -fbacktrace" 

CXX_FLAGSIC="-no-ipo -no-prec-div -heap-arrays -traceback" 
CXX_FLAGSPC="-traceback -Mdalign -gopt" 
CXX_FLAGSGC="-C " 

#
# We begin with null flags
#

FORT_FLAGS=""
CCFLAGS=""

#
# MKL link flags for intel, portland and gnu compilers
#

MKLIF_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl"
MKLPF_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl"
MKLGF_LIBS="-lmkl_gf_lp64    -lmkl_sequential -lmkl_core -ldl"

#
# No Default targeted architecture
#

ARCHIF=""
ARCHPF=""
ARCHGF=""


#
# special mkl module for fft not compiled by default
#

MKL_MOD=""
MKL_MOD_SOURCE=""

# default extension for binary files
EXEEXT=".prod"                                                    # Default is no extension (Automake variable)
#
# Enable debug mode
#
 AC_ARG_ENABLE([debug],
               [AS_HELP_STRING([--enable-debug], [Enable debug mode (check array bounds, implicit none, etc...). Should not be active in normal operations @<:@default is no@:>@])],
               [EXEEXT=".debug"                                                   # Extension is _debug (Automake variable)
                FCFLAGS="-O0 -g"                                                  # FCFLAGS is no optimization
                CCFLAGS="-O0 -g"                                                  # CCFLAGS is no optimization
                AC_FC_CHECK_BOUNDS([:],[AC_MSG_FAILURE([$CHECKBOUND_FAILMSG])])   # Modify FCFLAGS
                AC_MSG_RESULT([${GRN}$ac_cv_fc_check_bounds${CLR}])
                FORT_FLAGS="$ac_cv_fc_check_bounds"                               # Modify FORT_FLAGS
                case "x$REALFC" in
                   xifort    ) IMPFLAG="-u";             WARNFLAG="-warn all";
                               CIMPFLAG="-u";           CWARNFLAG="-warn" ;;
                   xpgfortran) IMPFLAG="-Mdclchk"        WARNFLAG="-Minform=inform";
                               CIMPFLAG="";             CWARNFLAG="-Minform=inform";;
                   xgfortran ) IMPFLAG="-fimplicit-none" WARNFLAG="-Wall";;
                esac
                AC_MSG_RESULT([Force implicit none with $GRN$IMPFLAG$CLR flag])
                FCFLAGS="$FCFLAGS $IMPFLAG"                                       # Modify FCFLAGS
                CCFLAGS="$CCFLAGS $CIMPFLAG"                                      # Modify CCFLAGS
                AC_MSG_RESULT([Turn on warnings with $GRN$WARNFLAG$CLR flag])
                FCFLAGS="$FCFLAGS $WARNFLAG"                                      # Modify FCFLAGS
                CCFLAGS="$CCFLAGS $CWARNFLAG"                                     # Modify CCFLAGS
                DEBUG_MODE="DEBUG"
               ],
               [FCFLAGS="-O3" ; CCFLAGS="-O3"                                     #  maximal optimization
                DEBUG_MODE="NORMAL"
               ]
          )
#
# Enable AVX512 optimization
#

#
# Skylake   architecture
#

AC_ARG_ENABLE([skylake],
              [AS_HELP_STRING([--enable-skylake], [Enable AVX512 Optimization for Skylake Processors @<:@default is no@:>@])],
              [ARCHIF="-xCORE-AVX512 -qopt-zmm-usage=high";
               ARCHGF="-march=skylake-avx512";
               ARCHPF="-tp=skylake"
              ]
             )
#
# KNL Xeon Phi architecture
#

AC_ARG_ENABLE([knl],
              [AS_HELP_STRING([--enable-knl], [Enable AVX512 Optimization for KNL (Xeon Phi) Processors @<:@default is no@:>@])],
              [ARCHIF="-xCOMMON-AVX512";
               ARCHGF="-march=knl";
               ARCHPF="-tp=knl"
              ]
             )

#
# Enable FFT modes (default is generic) FFTW_TYPE will contain the value chosen by the user
#

FFTW_TYPE=generic

AC_ARG_ENABLE([fft-generic],
              [AS_HELP_STRING([--enable-fft-generic], [Enable generic FFT mode @<:@default is yes@:>@])],
              [FFTW_TYPE=generic],
              [:]
             )
AC_ARG_ENABLE([fft-mkl],
              [AS_HELP_STRING([--enable-fft-mkl], [Enable MKL FFT mode @<:@default is no@:>@])],
              [FFTW_TYPE=mkl; with_fftlib=mkl],
              [:]
             )
AC_ARG_ENABLE([fft-fftw3],
              [AS_HELP_STRING([--enable-fft-fftw3], [Enable fftw3 FFT mode @<:@default is no@:>@])],
              [FFTW_TYPE=fftw3],
              [:]
             )
AC_ARG_ENABLE([fft-fftw3_f03],
              [AS_HELP_STRING([--enable-fft-fftw3_f03], [Enable fftw3_f03  FFT mode @<:@default is no@:>@])],
              [FFTW_TYPE=fftw3_f03],
              [:]
             )
#
# Enable plumed interface. PLUMED will contain the value chosen by the user
#

PLUMED=
PLUMEDFLAG=
LIBPLUMED_PATH=
LIBPLUMED_LIB=
LIBPLUMED_INC=
TCL=
TCLFLAG=
LIBTCL_PATH=
LIBTCL_LIB=
LIBTCL_INC=

AC_ARG_ENABLE([plumed],
              [AS_HELP_STRING([--enable-plumed], [Enable plumed interface @<:@default is no@:>@])],
              [PLUMED="PLUMED"],
              [:]
             )


#
# Enable colvars interface. COLVARS will contain the value chosen by the user
#

COLVARS=
COLVARSFLAG=
LIBCOLVARS_PATH=
LIBCOLVARS_LIB=
LIBCOLVARS_INC=

AC_ARG_ENABLE([colvars],
              [AS_HELP_STRING([--enable-colvars], [Enable colvars interface @<:@default is no@:>@])],
              [COLVARS="COLVARS"; TCL="TCL"],
              [:]
             )

#
# Enable tcl interface. TCL will contain the value chosen by the user
#


AC_ARG_ENABLE([tcl],
              [AS_HELP_STRING([--enable-tcl], [Enable tcl library @<:@default is no@:>@])],
              [TCL="TCL"],
              [:]
             )



#
# BLAS checking. BLAS will be MKL or LAPACK 
#

AC_MSG_CHECKING([what type of BLAS the user wants])
AC_ARG_WITH([blaslib],
            [AS_HELP_STRING([--with-blaslib=<BLAS LIB>],
                            [Specify BLAS library @<:@mkl, lapack or /absolute/path/to/BLAS_library@:>@])],
            [BLAS_LIB=$with_blaslib]
           )
AC_MSG_RESULT([${GRN}$BLAS_LIB${CLR}])

if test -z $BLAS_LIB; then BLAS_LIB=mkl;fi # If no BLAS type given, default to MKL

BLASUNKMSG="${RED}Unknown BLAS lib ${YEL}$BLAS_LIB${RED} : must be${GRN} mkl${RED} ,${GRN} lapack${RED}"
BLASUNKMSG="$BLASUNKMSG or a ${GRN}valid path to a BLAS library. ${CLR}"
case "$BLAS_LIB" in
    mkl|/*mkl*)                                               # mkl case (either mkl, or a path containing mkl)
       if test -z "$MKLROOT" ; then                           # MKLROOT should contains the root path for MKL libraries
          MKLROOT=$BLAS_LIB                                   # Could be defined by modules or user before configure begins
       fi
       AS_IF([test $(basename $MKLROOT) = intel64],
             [MKLROOT=$(dirname $(dirname $MKLROOT))],        # Remove the final /lib/intel64 to get the root directory of MKL
             [test $(basename $MKLROOT) = lib],
             [MKLROOT=$(dirname $MKLROOT)],                   # Remove the final /lib to get the root directory of MKL
             [:]
            )
       MKL_LIB_PATH="$MKLROOT/lib/intel64"
       MKL_INCLUDE_PATH="$MKLROOT/include"

       MKLIF="$MKLIF_LIBS"
       MKLPF="$MKLPF_LIBS"
       MKLGF="-Wl,--no-as-needed $MKLGF_LIBS"

       LDFLAGS="-L${MKL_LIB_PATH}";                            # Define the library path to be used by AC_CHECK_LIB
#
# Define the right link flags depending on the compiler
# These flags will be used by the AC_CHECK_LIB Macro
#
       case "x$REALFC" in
          xifort    ) MKL="$MKLIF" ; MKL_ADD_LIBS="$MKLIF_LIBS";;
          xpgfortran) MKL="$MKLPF" ; MKL_ADD_LIBS="$MKLPF_LIBS";;
          xgfortran ) MKL="$MKLGF" ; MKL_ADD_LIBS="$MKLGF_LIBS";;
       esac
       BLASFAILMSG="${RED}Cannot use $BLAS_LIB : missing functions or library ${YEL}mkl_core ${RED}not found${CLR}"

# 
# If the test is true, BLAS will contain the full link command
# to be exported in the different Makefiles
#
       AC_CHECK_LIB( [mkl_core],
                     [dgetrs],
                     [BLAS="$LDFLAGS $MKL"],
                     [AC_MSG_FAILURE([$BLASFAILMSG])],
                     [$MKL_ADD_LIBS]
                   )
       BLASTYPE="MKL"                                          # These 3 variables are just for the NOTICE macro
       BLAS_LIBS="$MKL"
       BLAS_LIBS_PATH="$MKL_LIB_PATH"
       ;;
    lapack|/*lapack*)                                          # lapack case (either lapack, or a path containing lapack)
       if test -z "$LAPACK" ; then                             # LAPACK should contains the root path for LAPACK libraries
          LAPACK=$BLAS_LIB                                     # Could be defined by modules or user before configure begins
       fi
      AS_IF([test $(basename $LAPACK) = lib],
            [LAPACK=$(dirname $LAPACK)],                       # Remove the final /lib to get the root directory of LAPACK
            [:]
           )
       LAPACK_LIB_PATH=$LAPACK
       LAPACKGF="-llapack"
       LAPACKIF="-llapack"
       LAPACKPF="-llapack"

       LAPACK_ADD_LIBS="-lrefblas"
       BLASFAILMSG="${RED}Cannot use $BLAS_LIB : missing functions or library ${YEL}lapack ${RED}not found${CLR}"
       LDFLAGS="-L ${LAPACK_LIB_PATH}";                        # Define the library path to be used by AC_CHECK_LIB
#
# Define the right link flags depending on the compiler
# These flags wiil be used by the AC_CHECK_LIB Macro
#
       case "x$REALFC" in
          xifort    ) LAPACK="$LAPACKIF $LAPACK_ADD_LIBS";; 
          xpgfortran) LAPACK="$LAPACKPF $LAPACK_ADD_LIBS";;
          xgfortran ) LAPACK="$LAPACKGF $LAPACK_ADD_LIBS";;
       esac
# 
# If the test is true, BLAS will contain the full link command
# to be exported in the different Makefiles
#
       AC_CHECK_LIB( [lapack],
                     [dgetrs],
                     [BLAS="$LDFLAGS $LAPACK"],
                     [AC_MSG_FAILURE([$BLASFAILMSG])],
                     [$LAPACK_ADD_LIBS]
                   )
       BLASTYPE="LAPACK"                                       # These 3 variables are just for the NOTICE macro
       BLAS_LIBS="$LAPACK"
       BLAS_LIBS_PATH="$LAPACK_LIB_PATH"
       ;;
   *)  AC_MSG_FAILURE( [$BLASUNKMSG])                          # No MKL, no lapack, die !!!
       ;;
esac

#
# FFT Checking
#

AC_MSG_CHECKING([what FFT Library the user wants])
AC_ARG_WITH([fftlib],
            [AS_HELP_STRING([--with-fftlib=<FFT LIB>],
                            [Specify a library for FFT called by 2decomp @<:@mkl or fftw3@ or /absolute/path/to/FFTW3_library:>@])],
            [FFTW_LIB=$with_fftlib]
           )
AC_MSG_RESULT([${GRN}$FFTW_LIB${CLR}])

if test -z $FFTW_LIB; then FFTW_LIB="none";fi                  # No FFTW_LIB given, so none

FFTUNKMSG="${RED}Unknown FFT LIB ${YEL}$FFTW_LIB${RED} : must be${GRN} mkl${RED} or${GRN} fftw3${RED} ${RED}.${CLR}"
case "$FFTW_LIB" in
    none )                                                     # none case
      AS_IF([test "x$FFTW_TYPE" = "xgeneric"],                 # This is legal, as long as we use the generic implementation in 2decomp 
            [FFTW_INCLUDES="-I ."; FFTW_LIB="" FFTW_INCLUDES_PATH="." ],
            [AC_MSG_FAILURE(["$NO_FFT_FAILMSG"])]              # If not, we absolutely need a fftw3 library
           )
       ;;
    mkl | /*mkl*)                                              # mkl case. Almost the same code as for BLAS
       if test -z "$MKLROOT" ; then
          MKLROOT=$FFTW_LIB
       fi
       AS_IF([test $(basename $MKLROOT) = intel64],
             [MKLROOT=$(dirname $(dirname $MKLROOT))],         # Remove the final /lib/intel64 to get the root directory of MKL
             [test $(basename $MKLROOT) = lib],
             [MKLROOT=$(dirname $MKLROOT)],                    # Remove the final /lib to get the root directory of MKL
             [:]
            )
       MKL_LIB_PATH="$MKLROOT/lib/intel64"
       MKL_INCLUDE_PATH="$MKLROOT/include"
       FFTW_INCLUDES_PATH="$MKL_INCLUDE_PATH"
       FFTW_MKL_INCLUDES_PATH="$MKL_INCLUDE_PATH"
       MKLIF=" $MKLIF_LIBS"
       MKLPF=" $MKLPF_LIBS"
       MKLGF="-Wl,--no-as-needed $MKLGF_LIBS"
#
# Define the right link flags depending on the compiler
# These flags wiil be used by the AC_CHECK_LIB Macro
#
       case "x$REALFC" in
          xifort    ) MKL="$MKLIF" ; MKL_ADD_LIBS="$MKLIF_LIBS";;
          xpgfortran) MKL="$MKLPF" ; MKL_ADD_LIBS="$MKLPF_LIBS";;
          xgfortran ) MKL="$MKLGF" ; MKL_ADD_LIBS="$MKLGF_LIBS";;
       esac
       FFTFAILMSG="${RED}Cannot use $FFTW_TYPE : missing functions or library ${YEL}mkl_core ${RED}not found${CLR}"
       LDFLAGS="-L${MKL_LIB_PATH}";
       AC_CHECK_LIB([mkl_core],[dfftw_plan_dft_c2r_2d],         #  We test a specific fftw function of MKL
                     [:],
                     [AC_MSG_FAILURE([$FFTFAILMSG])],
                     [$MKL_ADD_LIBS])
       FFTW_INCLUDES="-I $MKL_INCLUDE_PATH"
       FFTW_LIBS="$LDFLAGS$MKL"
       FFTWLIB="$MKL"
       FFTW_LIB_PATH="$MKL_LIB_PATH"
       AS_IF([test "x$FFTW_TYPE" = "xfftw3" || test "x$FFTW_TYPE" = "xfftw3_f03"], # If we want the fftw3* interface of MKL
              [FFTW_INCLUDES_PATH="${FFTW_INCLUDES_PATH}/fftw";                    # Add fftw to the include path
               FFTW_INCLUDES="${FFTW_INCLUDES}/fftw"],
              [ MKL_MOD_BIN="mkl_dfti.o"                                           # if not, compile the special mkl module
                MKL_MOD_SOURCE="${MKL_MOD_BIN%%.o}.f90"
              ]
            )
       FFTWTYPE="MKL"                                           # This is just for the NOTICE
       ;;
    fftw3 | /*fftw*)                                            # fftw3 case
       if test -z "$FFTW" ; then
          FFTW=$FFTW_LIB
       fi
       AS_IF([test $(basename $FFTW) = lib],
             [FFTW=$(dirname $FFTW)],                           # Remove the final /lib to get the root directory of FFTW3
             [:]
            )
       FFTW_LIB_PATH="$FFTW/lib"
       FFTW_INCLUDES_PATH="$FFTW/include"
       FFTWIF="-lfftw3"
       FFTWPF="-lfftw3"
       FFTWGF="-lfftw3"
       FFTW_ADD_LIBS=""
       FFTWFAILMSG="${RED}Cannot use $FFTW_TYPE : missing functions or library ${YEL} fftw3 ${RED}not found${CLR}"
       LDFLAGS="-L ${FFTW_LIB_PATH}";
#
# Define the right link flags depending on the compiler
# These flags wiil be used by the AC_CHECK_LIB Macro
#
# Define also the right preprocessor flags depending on the compiler
#
       case "x$REALFC" in
          xifort    ) FFTWLIB="$FFTWIF";;
          xpgfortran) FFTWLIB="$FFTWPF";;
          xgfortran ) FFTWLIB="$FFTWGF";;
       esac
#
# Test the specific fftw3 function
#

       AC_CHECK_LIB([fftw3],[dfftw_plan_dft_c2r_2d],
                    [FFTW_INCLUDES="-I $FFTW_INCLUDES_PATH"
                     FFTW_LIBS="$LDFLAGS $FFTWLIB"],
                    [AC_MSG_FAILURE([$FFTWFAILMSG])],
                    [$FFTW_ADD_LIBS])
       FFTWTYPE="FFTW3"                                         # This is just for the NOTICE
       ;;
    *)  AC_MSG_FAILURE( [$FFTUNKMSG]) ;;
esac

#
# setup plumed 
#
PLUMBADCOMPILMSG="${RED} Invalid wrapped compiler ${YEL}$REALFC${RED} for PLUMED : must be ${GRN} ifort ${RED}.${CLR}"

       AS_IF([test x$PLUMED = xPLUMED],
             [ PLUMEDFLAG="-DPLUMED"
               case "x$REALFC" in
                  xifort    ) COMP="Intel"; CXXCOMP="mpiicpc";;
                  xpgfortran) COMP="Portland"; CXXCOMP="mpic++";;
                  xgfortran ) COMP="Gnu"; CXXCOMP="mpic++";;
               esac
               LIBPLUMED_PATH="$PREFIX/plumed/$COMP"
               LIBPLUMED_LIBPATH="$LIBPLUMED_PATH/lib/"
               LIBPLUMED_LIB="-L $LIBPLUMED_LIBPATH -lplumed -lplumedKernel -lstdc++"
               LIBPLUMED_INC="-I $LIBPLUMED_PATH/include"
               PLUEXT="_plumed"
               PLUMEDINCLUDE="plumed.mk"
               PLUMEDDEP="$LIBPLUMED_LIBPATH/libplumed.so $LIBPLUMED_LIBPATH/libplumedKernel.so"
               AC_MSG_NOTICE(["$YEL Configuring PLUMED in $LIBPLUMED_PATH$GRN${CLR}"])
               #
               # States that this Makefiles.conf should be generated
               #
               AC_CONFIG_FILES(plumed2/Makefile.conf, [(cd plumed2 ; ./configure MPICXX=$CXXCOMP --prefix=$PREFIX  --enable-mpi)],[PREFIX=$(pwd)/plumed/$COMP/; CXXCOMP=$CXXCOMP ] )
             ],
             [ PLUMEDFLAG=""
               LIBPLUMED_PATH=""
               LIBPLUMED_LIB=""
               LIBPLUMED_INC="" 
               PLUEXT=""
               PLUMEDINCLUDE="noplumed.mk"
               PLUMEDDEP=""
	     ]
            )
#
# setup colvars 
#
COLVBADCOMPILMSG="${RED} Invalid wrapped compiler ${YEL}$REALFC${RED} for COLVARS : must be ${GRN} ifort ${RED}.${CLR}"
       AS_IF([test x$COLVARS = xCOLVARS],
             [ COLVARSFLAG="-DCOLVARS -DCOLVARS_TCL"
               case "x$REALFC" in
                  xifort    ) COMP="Intel";;
                  xpgfortran) COMP="Portland";;
                  xgfortran ) COMP="Gnu";;
               esac
               LIBCOLVARS_PATH="../colvars/"
               LIBCOLVARS_LIBPATH="$LIBCOLVARS_PATH/"
               LIBCOLVARS_LIB="-L $LIBCOLVARS_LIBPATH/ -lcolvars -lstdc++ "
               LIBCOLVARS_INC="-I $LIBCOLVARS_PATH -I $PREFIX/tcl8.16.3/generic"
               LIBCOLVARS_CPPFLAGS="$LIBCOLVAR_INC"
               COLVEXT="_colvars"
               COLVARSINCLUDE="colvars.mk"
               TINKERCOLVARS="libtinkerclv_a"
               COLVARSDEP="$LIBCOLVARS_PATH/libcolvars.a"
             ],
             [ COLVARSFLAG=""
               LIBCOLVARS_PATH=""
               LIBCOLVARS_LIB=""
               LIBCOLVARS_INC="" 
               COLVEXT=""
               COLVARSINCLUDE="nocolvars.mk"
               TINKERCOLVARS=""
               COLVARSDEP=""
	     ]
            )
#
# setup tcl 
#
TCLBADCOMPILMSG="${RED} Invalid wrapped compiler ${YEL}$REALFC${RED} for TCL : must be ${GRN} icc ${RED}.${CLR}"

       AS_IF([test x$TCL = xTCL],
             [ TCLFLAG="-DTCL"
               case "x$REALFC" in
                  xifort    ) COMP="Intel"; CCOMP="icc";;
                  xpgfortran) COMP="Portland"; CCOMP="pgcc";;
                  xgfortran ) COMP="Gnu"; CCOMP="gcc";;
               esac
               LIBTCL_EXEPATH="$PREFIX/tcl8/$COMP"
               LIBTCL_EXELIBPATH="$LIBTCL_EXEPATH/lib/"
               LIBTCL_PATH="$PREFIX/tcl8.6.13/unix/"
               LIBTCL_LIBPATH="$LIBTCL_PATH"
               LIBTCL_LIB="-L $LIBTCL_LIBPATH -ltcl8.6"
               LIBTCL_INC="-I $LIBTCL_PATH/../generic"
               TCLEXT="_tcl"
               TCLINCLUDE="tcl.mk"
               AC_MSG_NOTICE(["$YEL Configuring TCL in $LIBTCL_PATH$GRN${CLR}"])
               #
               # States that this Makefiles.conf should be generated
               #
               AC_CONFIG_FILES(tcl8.6.13/unix/Makefile.conf, [(cd $LIBTCL_PATH ; ./configure CC=$CCOMP --prefix=$PREFIXTCL)],[PREFIXTCL=$(pwd)/tcl8/$COMP/ ; CCOMP=$CCOMP ; LIBTCL_PATH=$LIBTCL_PATH] )
             ],
             [ TCLFLAG=""
               LIBTCL_PATH=""
               LIBTCL_LIB=""
               LIBTCL_INC="" 
               TCLEXT=""
               TCLINCLUDE="notcl.mk"
	     ]
            )
#
# Final setting of the Compiler flags
#
       case "x$REALFC" in
          xifort    ) FORT_FLAGS="$FORT_FLAGS $FORT_FLAGSIF $ARCHIF";PREPROFLAG="-fpp";
                      CXX_FLAGS="$CXX_FLAGSIC $ARCHIF";;
          xpgfortran) FORT_FLAGS="$FORT_FLAGS $FORT_FLAGSPF $ARCHPF";PREPROFLAG="-Mpreprocess";
                      NVMPICXX="-lmpi_cxx"; # Added in LDADD when using pgfortran
                      CXX_FLAGS="$CXX_FLAGSPC $ARCHPF";;
          xgfortran ) FORT_FLAGS="-Wno-argument-mismatch $FORT_FLAGS $FORT_FLAGSGF $ARCHGF";PREPROFLAG="-cpp";
                      CXX_FLAGS="$CXX_FLAGSGC $ARCHGF";;
       esac

#
# Flags for .cpp sources 
#
CXXFLAGS="-std=c++11 $COLVARSFLAG $PLUMEDFLAG $TCLFLAG $CCFLAGS $CXX_FLAGS"
AM_CPPFLAGS="$LIBCOLVARS_CPPFLAGS"

#
# Flags for .f sources are the same as those for .f90 sources
# Fortran compiler is the wrapper
FCFLAGS="$PREPROFLAG $COLVARSFLAG $PLUMEDFLAG $TCLFLAG $FCFLAGS $FORT_FLAGS"
FFLAGS=$FCFLAGS
FC=$MPIFC

#
# if we have MKL for both BLAS and FFTW3, just remove the FFTW3 stuff from the link command to avoid duplication
#

AS_IF([test "x$FFTW_LIBS" = "x$BLAS" ],
      [FFTW_LIBS=]
     )
#
# Set the LDFLAGS to nothing. We don't need it in the Makefiles, because we built all the link commands before
#

LDFLAGS=

#
# Build the messages for DEBUG mode
#

AS_IF([test "x$DEBUG_MODE" = xDEBUG ],
      [RUNMO_MSG="** Running Mode         : ${YEL}$PLUMED $COLVARS $TCL $DEBUG_MODE$GRN (binaries'extension is $YEL'$EXEEXT'$GRN)"],
      [RUNMO_MSG="** Running Mode         : ${YEL}$PLUMED $COLVARS $TCL $DEBUG_MODE$GRN (No extension for binaries)"]
     )
#
# set the extension of binaries
# composite of _plumed, _colvars and .debug or .prod 
#
EXEEXT="$PLUEXT$COLVEXT$TCLEXT$EXEEXT"

#
# Build the messages for NOTICE
#

HEAD1_MSG="**********************************************************************"
HEAD2_MSG="**"

#RUNMO_MSG="** ${CLR}Running Mode         : $DEBUG_MSG"
MPICO_MSG="** MPI Fortran Wrapper  : ${GRN}$MPIFC"
FORCO_MSG="** Fortran Compiler     : ${GRN}$REALFC"
FOFL1_MSG="** Fortran flags        : ${GRN}${FCFLAGS:0:73}"
FOFL2_MSG="**                      : ${GRN}${FCFLAGS:73:146}"
DECLI_MSG="** 2decomp Library      : ${GRN}$LIB2DECOMP_LIB"
PLULI_MSG="** PLUMED  Library      : ${GRN}$LIBPLUMED_LIB"
PLUIN_MSG="** PLUMED  Includes     : ${GRN}$LIBPLUMED_INC"
COLLI_MSG="** COLVARS  Library     : ${GRN}$LIBCOLVARS_LIB"
COLIN_MSG="** COLVARS  Includes    : ${GRN}$LIBCOLVARS_INC"
TCLLI_MSG="** TCL  Library         : ${GRN}$LIBTCL_LIB"
TCLIN_MSG="** TCL  Includes        : ${GRN}$LIBTCL_INC"
FFTTY_MSG="** FFTW3 Interface      : ${GRN}$FFTW_TYPE of the $FFTWTYPE library"
FFTPA_MSG="** FFTW3 Path           : ${GRN}$FFTW_LIB_PATH"
FFTIN_MSG="** FFTW3 Includes       : ${GRN}$FFTW_INCLUDES"
FFTLI_MSG="** FFTW3 Library        : ${GRN}$FFTWLIB"
BLATY_MSG="** BLAS Type            : ${GRN}$BLASTYPE"
BLAPA_MSG="** BLAS Path            : ${GRN}$BLAS_LIBS_PATH"
PREFI_MSG="** Prefix installation  : ${GRN}$PREFIX"
BINAR_MSG="** Binaries location    : ${GRN}$PREFIX/bin"

RUNMO_MSG="$RUNMO_MSG${PADD:0:$((120 - ${#RUNMO_MSG}))}$CLR"
MPICO_MSG="$CLR$MPICO_MSG${PADD:0:$((120 - ${#MPICO_MSG}))}$CLR"
FORCO_MSG="$CLR$FORCO_MSG${PADD:0:$((120 - ${#FORCO_MSG}))}$CLR"
FOFL1_MSG="$CLR$FOFL1_MSG${PADD:0:$((120 - ${#FOFL1_MSG}))}$CLR"
FOFL2_MSG="$CLR$FOFL2_MSG${PADD:0:$((120 - ${#FOFL2_MSG}))}$CLR"
DECLI_MSG="$CLR$DECLI_MSG${PADD:0:$((120 - ${#DECLI_MSG}))}$CLR"
PLULI_MSG="$CLR$PLULI_MSG${PADD:0:$((120 - ${#PLULI_MSG}))}$CLR"
PLUIN_MSG="$CLR$PLUIN_MSG${PADD:0:$((120 - ${#PLUIN_MSG}))}$CLR"
COLLI_MSG="$CLR$COLLI_MSG${PADD:0:$((120 - ${#COLLI_MSG}))}$CLR"
COLIN_MSG="$CLR$COLIN_MSG${PADD:0:$((120 - ${#COLIN_MSG}))}$CLR"
TCLLI_MSG="$CLR$TCLLI_MSG${PADD:0:$((120 - ${#TCLLI_MSG}))}$CLR"
TCLIN_MSG="$CLR$TCLIN_MSG${PADD:0:$((120 - ${#TCLIN_MSG}))}$CLR"
FFTTY_MSG="$CLR$FFTTY_MSG${PADD:0:$((120 - ${#FFTTY_MSG}))}$CLR"
FFTPA_MSG="$CLR$FFTPA_MSG${PADD:0:$((120 - ${#FFTPA_MSG}))}$CLR"
FFTIN_MSG="$CLR$FFTIN_MSG${PADD:0:$((120 - ${#FFTIN_MSG}))}$CLR"
FFTLI_MSG="$CLR$FFTLI_MSG${PADD:0:$((120 - ${#FFTLI_MSG}))}$CLR"
BLATY_MSG="$CLR$BLATY_MSG${PADD:0:$((120 - ${#BLATY_MSG}))}$CLR"
BLAPA_MSG="$CLR$BLAPA_MSG${PADD:0:$((120 - ${#BLAPA_MSG}))}$CLR"
PREFI_MSG="$CLR$PREFI_MSG${PADD:0:$((120 - ${#PREFI_MSG}))}$CLR"
BINAR_MSG="$CLR$BINAR_MSG${PADD:0:$((120 - ${#BINAR_MSG}))}$CLR"

#
# Make the substitutions in all Makefiles.am, Makefiles.in and Makefiles.
#

AC_SUBST([FCFLAGS])
AC_SUBST([CPPFLAGS])
AC_SUBST([MKL_MOD_BIN])
AC_SUBST([MKL_MOD_SOURCE])
AC_SUBST([LIB2DECOMP_INC])
AC_SUBST([LIB2DECOMP_LIB])
AC_SUBST([PLUMED])
AC_SUBST([LIBPLUMED_INC])
AC_SUBST([LIBCOLVARS_INC])
AC_SUBST([LIBPLUMED_LIB])
AC_SUBST([LIBPLUMED_LIBPATH])
AC_SUBST([PLUMEDDEP],[$PLUMEDDEP])
AC_SUBST([TCL])
AC_SUBST([LIBTCL_INC])
AC_SUBST([LIBTCL_LIB])
AC_SUBST([LIBTCL_LIBPATH])
AC_SUBST([LIBTCL_EXELIBPATH])
AC_SUBST([FFTW_TYPE])
AC_SUBST([FFTW_LIB_PATH])
AC_SUBST([FFTW_INCLUDES])
AC_SUBST([FFTW_INCLUDES_PATH])
AC_SUBST([FFTW_MKL_INCLUDES_PATH])
AC_SUBST([FFTW_LIBS])
AC_SUBST([BLAS])
AC_SUBST([NVMPICXX])
AC_SUBST([PYTHON])
AC_SUBST([PREFIX])
AC_SUBST([INCLUDE_DEPEND_MK],["include depend.mk"])
AC_SUBST([INCLUDE_PLUMED_MK],["include "$PLUMEDINCLUDE])
AC_SUBST([INCLUDE_TCL_MK],["include "$TCLINCLUDE])
AC_SUBST([INCLUDE_COLVARS_MK],["include "$COLVARSINCLUDE])
AC_SUBST([LIBTINKERCLV],[$TINKERCOLVARS])
AC_SUBST([LIBCOLVARS_LIB],[$LIBCOLVARS_LIB])
AC_SUBST([COLVARSDEP],[$COLVARSDEP])

#
# States that these Makefiles should be generated
#

AC_CONFIG_FILES( Makefile
                 2decomp_fft/src/Makefile
                 example/Makefile
                 params/Makefile
                 source/Makefile
                 colvars/Makefile
                 tutorials/Makefile
                 tutorials/SMD/Makefile
                 tutorials/SMD/CFSMD/Makefile
                 tutorials/SMD/CVSMD/Makefile
                 scripts/Makefile
		 scripts/tinkerhp
               )

#
# Produces  Makefiles.in and Makefiles.
#

AC_OUTPUT

#
# Fancy (well, I think...) output for the user
#
AC_MSG_NOTICE()
AC_MSG_NOTICE([$HEAD1_MSG])
AC_MSG_NOTICE([$HEAD2_MSG])
AC_MSG_NOTICE([$RUNMO_MSG])
AC_MSG_NOTICE([$MPICO_MSG])
AC_MSG_NOTICE([$FORCO_MSG])
AC_MSG_NOTICE([$FOFL1_MSG])
AC_MSG_NOTICE([$FOFL2_MSG])
AC_MSG_NOTICE([$DECLI_MSG])
AC_MSG_NOTICE([$PLULI_MSG])
AC_MSG_NOTICE([$PLUIN_MSG])
AC_MSG_NOTICE([$COLLI_MSG])
AC_MSG_NOTICE([$COLIN_MSG])
AC_MSG_NOTICE([$TCLLI_MSG])
AC_MSG_NOTICE([$TCLIN_MSG])
AC_MSG_NOTICE([$FFTTY_MSG])
AC_MSG_NOTICE([$FFTPA_MSG])
AC_MSG_NOTICE([$FFTIN_MSG])
AC_MSG_NOTICE([$FFTLI_MSG])
AC_MSG_NOTICE([$BLATY_MSG])
AC_MSG_NOTICE([$BLAPA_MSG])
AC_MSG_NOTICE([$PREFI_MSG])
AC_MSG_NOTICE([$BINAR_MSG])
AC_MSG_NOTICE([$HEAD2_MSG])
AC_MSG_NOTICE([$HEAD1_MSG])
AC_MSG_NOTICE()
