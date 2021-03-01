#!/bin/bash

# This script install proceed to insatallation of Tinker-HP
# Some help on makefile's target can be found inside build.md in the root directory

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
   !!! WARNING GNUROOT variable not find in environment !!!
   setting it to "$GNUROOT" according to g++ compiler location

   *** You might later need to export this variable for development purposes ***

   Installation will resume shortly
EOP
   sleep 5
fi

#test if sourced or not
(return 0 2>/dev/null) && sourced=1 || sourced=0
if [ $sourced -eq 1 ]; then
   cat << END

ci/install : This script is not destined to be sourced
             Please Type \`bash ci/install.sh\` or \`ci/install.sh\`

END

else

#set -x
tinkerdir=$(dirname `dirname $0`)

# Carefully Add 'var=value' in common_config to customize your build
# Make sure those params do not conflict with 'common_config_[dm]'
c_c=60,70
cuda_ver=10.1
build_plumed=0

common_config="compute_capability=$c_c cuda_version=$cuda_ver PLUMED_SUPPORT=$build_plumed"
common_config_d="$common_config prog_suffix=.gpu"
common_config_m="$common_config prec=m prog_suffix=.mixed"

[ $# -ge 1 ]&& par=$1 || par=6

cd $tinkerdir && \
mkdir -p bin && \
ln -sf Makefile.pgi source/Makefile  && \
cd source

cat << END

------ Compiling 2decomp_fft library ------

END
make $common_config_d 2decomp_fft

cat << END

------ Compiling thrust wrapper      ------

END
make $common_config thrust

if [ $build_plumed -eq 1 ]; then
   cat << END

------ Building PLUMED ------

END
   make plumed -j$par
fi

cat << END

====== Compiling TINKER-HP           ======

END
make $common_config_d -j$par
if [ "$?" != "0" ]; then
   cat << END

   CAUTION! Something went wrong during compilation !!"
            Please Fix the issue and run ci/install.sh again"

END
   cd ../
   exit
fi # Compiling test

echo '***  Cleaning objects files and modules ***'
make clean >/dev/null && \
cat << END

------ Compiling 2decomp_fft library in single precision ------

END
make $common_config_m 2decomp_fft_rebuild_single

cat << END

===== Recompiling TINKER-HP in mixed precision =====

END
make $common_config_m -j$par
if [ "$?" != "0" ]; then
   cat << END
             ------ WARNING ------
   Something went wrong during mixed precision compilation !!
   Please Fix the issue and enter following commands to resume
$> cd source
$> make $common_config_m -j$par
             ---------------------
END
   cd ../
   exit
fi # Mixed Compiling test

# Finalize
make clean >/dev/null && \
cat << END

  *****  TINKER-HP GPU installation is successfully completed  *****

END

fi # Source test
