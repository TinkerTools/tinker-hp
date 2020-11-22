#!/bin/bash

# This script install proceed to insatallation of Tinker-hP on dgx2 server
# Some help on makefile's target are written at the end of this file

# Make sure you have MKLROOT and GNUROOT variables in your environnement path before running

if [[ -z ${MKLROOT+x} ]]; then
   cat << EOP
   Installation ERROR !
   MKLROOT is unset in your environment !!!
   Please export That variable to MKL libraty Root and rerun this script

EOP
   exit
fi

if [[ -z ${GNUROOT+x} ]]; then
   GNUROOT=`which g++ | sed 's_/bin/g++__'`
   export GNUROOT
   cat << EOP
   !!! WARNING GNUROOT variable not find in env !!!
   setting it to "$GNUROOT" according to g++ compiler location

   *** You might later need to export this variable later in your environment ***

   Installation will resume shortly
EOP
   sleep 5
fi

#set -x
tinkerdir=`pwd`

[ $# -ge 1 ]&& par=$1 || par=4

ln -sf $tinkerdir/2decomp_fft/src/Makefile.inc.pgi $tinkerdir/2decomp_fft/src/Makefile.inc
ln -sf Makefile.pgi source/Makefile
cd source
cat << END

------ Compiling 2decomp_fft library ------

END
make 2decomp_fft_rebuild
cat << END

------ Compiling thrust wrapper      ------

END
make thrust_lib_rebuild
mkdir -p ../bin
cat << END

====== Compiling TINKER-HP           ======

END
make prog_suffix=.gpu -j$par
res=$?
#echo $res
if [ "$res" != "0" ]; then
   echo
   echo
   echo " CAUTION!  Something went wrong during compilation !!"
   echo "           Please Fix the issue and run ci/install.sh again"
   echo
   echo
   cd ../
else
make clean
cat << END

------ Compiling 2decomp_fft library in single precision ------

END
make 2decomp_fft_rebuild_single
cat << END

===== Recompiling TINKER-HP in mixed precision =====

END
make prec=m prog_suffix=.mixed -j$par
res=$?
if [ "$res" != "0" ]; then
   echo
   echo
   echo "Something went wrong during mixed precision compilation !!"
   echo "Please Fix the issue and enter following commands"
   echo "$> cd source"
   echo "$> make prec=m prog_suffix=.mixed -j6"
   echo
else
make clean
cat << END

TINKER-HP GPU intallation completes successully

   You may now proceed to simulations

END

fi
fi


###########################
# Some Help on Makefile use
###########################
: '
Some others targets available in _source/makefile_

  - create_build
    In association with _BUILD_DIR_, it is possible to make a new build directory
    with links to source files and _source/Makefile.pgi_. The default directory is `build/`
    e.g.  `make create_build BUIL_DIR=../build-mix`

  - 2decomp_fft_rebuild
    This phony target allow you to rebuild 2decomp_fft library in double precision
    Since the compilation uses the same directory for 2decomp_fft, you might want to
    be sure that your library is in the right precision

  - 2decomp_fft_rebuild_single
    Same as previous but in single precision

  - thrust_lib_rebuild
    Rebuild the wrapper on CUDA thrust library

  - all
    compile link and install `bin/dynamic` `bin/analyze` and `bin/minimize`

  - dynamic.mixed
    compile link and install </dynamic.mixed>

_______________________________

Makefile s configuration variables

  - arch  (host || device)
       allow you to select build target architecture '
