#!/bin/bash

# Tinker-HP [GPU] Quick Installation script
# Some help on makefile's target can be found inside build.md in the root directory

#test if sourced or not
(return 0 2>/dev/null) && sourced=1 || sourced=0
if [ $sourced -eq 1 ]; then
   cat << END

ci/install : This script is not destined to be sourced
             Please Type \`bash ci/install.sh\` or \`ci/install.sh\`

END

else

#set -x
tinkerdir=$(realpath $(dirname `dirname $0`))

#echo "Tinker-HP dir :" $tinkerdir
source $tinkerdir/ci/utils.sh
unset CC CXX FC

#
#=======================
# Configuration  Section

# Change those 'value' according to comments ahead to configure your build
#

dbuild=$tinkerdir/build0      #   Out-of-source build directory
target_arch='gpu'             #   [gpu]|cpu
c_c=60,70,80,86               #   Target GPU compute capability  [https://en.wikipedia.org/wiki/CUDA]
cuda_ver=11.7                 #   Cuda Version to be used by OpenACC compiler  (not the CUDA C/C++ compiler)
FPA=1                         #   Enable Fixed Precision Arithmetic (Useful for non HPC-Accelerators)
build_plumed=0                #   [0]|1     0: disable 1: enable   PLUMED  Interface
build_colvars=0               #   [0]|1     0: disable 1: enable   COLVARS Interface
NN=0                          #   [0]|1     0: disable 1: enable   Neural Network Python Interface
#add_host_f='-Mx,231,0x1'      #   Uncomment this when building Nvidia HPC-SDK package version 21.[3-7]

#  End Config
#=======================

current_config="arch=$target_arch compute_capability=$c_c cuda_version=$cuda_ver PLUMED_SUPPORT=$build_plumed COLVARS_SUPPORT=$build_colvars NN_SUPPORT=$NN"
[ -n "$add_options_f" ] && current_config="$current_config add_options_f=$add_options_f"
current_config_d="$current_config"
current_config_m="$current_config FPA_SUPPORT=$FPA prec=m"

# ------------------------------------
# Clean Project and exit if instructed
# ------------------------------------
[ $# -ge 1 ] && [ $1 = "clean" ] && cd $dbuild && make $current_config distclean && exit

[ $# -ge 1 ] && [ $1 -gt 1 ] && ntask=$1 || ntask=16
[ $# -ge 2 ] && [ $2 -gt 1 ] && ntask=$2 || ntask=16

chk_gnuroot

# ------------------------------------------------
# Check for MKLROOT variable in your environnement
# ------------------------------------------------
[[ -z ${MKLROOT+x} ]] && [[ ${target_arch} = 'cpu' ]] && error_mkl && exit

# --------
# Initiate
# --------
cd $tinkerdir && \
mkdir -p bin && \
[ ! -d $dbuild ] && cd source && \
make -f Makefile.pgi create_build BUILD_DIR=$dbuild
cd $dbuild

# ----------
# Build Step
# ----------
if [ $build_plumed -eq 1 ]; then
   in_notif "Building PLUMED" &&\
   make plumed -j$ntask
fi

# Double precision build
if [ $# -eq 0 ] || [ $# -eq 1 ] && [ "$1" != "nd" ]; then
   in_notif "Compiling TINKER-HP" && \
   make $current_config_d -j$ntask
   [ "$?" != "0" ] && error1st && exit # Compiling test

   in_notif 'Cleaning objects files and modules' && \
   make clean >/dev/null
fi

# Mixed precision build
[ ${FPA} -eq 1 ] && in_notif "Recompiling TINKER-HP in mixed precision + FPA support"
[ ${FPA} -eq 0 ] && in_notif "Recompiling TINKER-HP in mixed precision"
make $current_config_m 2decomp_fft_rebuild   && \
make $current_config_m -j$ntask dynamic bar pibar pimd
[ "$?" != "0" ] && error2nd  && exit # Mixed Compiling test

# --------
# Finalize
# --------
in_notif 'Cleaning object files and modules'         && \
make clean 2decomp_fft_clean thrust_clean >/dev/null && \

in_notif 'TINKER-HP GPU installation completes successfully'

fi # Source test
