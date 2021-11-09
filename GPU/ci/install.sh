#!/bin/bash

# Tinker-HP [GPU] Quick Installation script
# Some help on makefile's target can be found inside build.md in the root directory

# --------------------
# Function Definitions
# --------------------
in_notif(){
   printf "\n <<<<<  %-60s  >>>>> \n\n" "$1" && sleep 1
}
error1st(){
   cat << END

             ------ WARNING ------
   Something went wrong during compilation procedure "
   Please Fix the issue and run ci/install.sh again"
             ---------------------

END
   cd ..
}
error2nd(){
   cat << END

             ------ WARNING ------
   Something went wrong during mixed precision compilation !!
   Please Fix the issue and enter following commands to resume
$> cd source
$> make $current_config_m -j$ntask
             ---------------------

END
   cd ../
}
error_mkl(){
   cat << EOP

   !!! ci/install.sh  ERROR !!!
   >   MKLROOT is unset in your environment 
   --  Please export That variable to MKL library's Root and rerun this script

EOP
}


# Check for GNUROOT variable in your environnement path before running
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

# -------------
# Configuration
# -------------
#
# Carefully Add or Change 'var=value' in current_config to configure your installation
# Make sure those params do not conflict with 'current_config_[dm]'
#

c_c=60,70,80                  #   Target GPU compute capability  [https://en.wikipedia.org/wiki/CUDA]
cuda_ver=11.0                 #   Cuda Version to be used by OpenACC compiler  (not the CUDA C/C++ compiler)
build_plumed=0                #   [0]|1      0: disable 1: enable
target_arch='gpu'             #   [gpu]|cpu
FPA=1                         #   Enable Fixed Precision Arithmetic (Useful for non HPC-Accelerators)
#add_options_f='-Mx,231,0x1'   #   Uncomment this when building Nvidia HPC-SDK package strictly above 21.3 version

[ $FPA -eq 1 ] && p_suffix=fpa || p_suffix=mixed       # Config binary suffix following FPA feature
current_config="compute_capability=$c_c cuda_version=$cuda_ver PLUMED_SUPPORT=$build_plumed arch=$target_arch"
[ -n "$add_options_f" ] && current_config="$current_config add_options_f=$add_options_f"

current_config_d="$current_config prog_suffix=.gpu"
current_config_m="$current_config FPA_SUPPORT=$FPA prec=m prog_suffix=.$p_suffix"

[ $# -ge 1 ] && ntask=$1 || ntask=12

# Check for MKLROOT variable in your environnement
[[ -z ${MKLROOT+x} ]] && [[ ${target_arch} = 'cpu' ]] && error_mkl && exit

# --------
# Initiate
# --------
cd $tinkerdir && \
mkdir -p bin && \
ln -sf Makefile.pgi source/Makefile  && \
cd source

# -----
# Build
# -----
if [ $build_plumed -eq 1 ]; then
   in_notif "Building PLUMED" && make plumed -j$ntask
fi

in_notif "Compiling TINKER-HP" && \
make $current_config_d -j$ntask
[ "$?" != "0" ] && error1st && exit # Compiling test

in_notif 'Cleaning objects files and modules' && \
make clean >/dev/null

[ ${FPA} -eq 1 ] && in_notif "Recompiling TINKER-HP in mixed precision + FPA support"
[ ${FPA} -eq 0 ] && in_notif "Recompiling TINKER-HP in mixed precision"
make $current_config_m 2decomp_fft_rebuild_single   && \
make $current_config_m -j$ntask
[ "$?" != "0" ] && error2nd  && exit # Mixed Compiling test

# --------
# Finalize
# --------
in_notif 'Cleaning object files and modules'         && \
make clean 2decomp_fft_clean thrust_clean >/dev/null && \

in_notif 'TINKER-HP GPU installation completes successfully'

fi # Source test
