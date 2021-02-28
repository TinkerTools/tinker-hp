if [ -n ${TINKERHP_ROOT} ]; then
   export PATH="$TINKERHP_ROOT/plumed/lib/:$PATH"
   export LIBRARY_PATH="$TINKERHP_ROOT/plumed/lib/:$LIBRARY_PATH"
   export LD_LIBRARY_PATH="$TINKERHP_ROOT/plumed/lib/:$LD_LIBRARY_PATH"
   export PLUMED_KERNEL="$TINKERHP_ROOT/plumed/lib/libplumedKernel.so"
   export PLUMED_VIMPATH="$TINKERHP_ROOT/plumed2/vim"
   export PYTHONPATH="$TINKERHP_ROOT/plumed2/python:$PYTHONPATH"
else
   echo '             **********************************                 '
   echo ' ***************************************************************'
   echo '                         WARNING !!!   '
   echo ' {TINKERHP_ROOT} is missing or empty from your environment !'
   echo
   echo ' Please provide that variable to point Tinker-HP main directory '
   echo ' before sourcing this file '
   echo ' ***************************************************************'
   echo '             **********************************                 '
fi
