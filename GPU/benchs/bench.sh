#!/bin/bash

source ../scripts.sh  # Load functions

#exe=../../bin/dynamic         # Double precision executable
exe=../../bin/dynamic.mixed   # Mixed Precision executable (Single float| Double Precision)

[ -n "$1" ] && system=$1
[ -z "$system" ] && echo 'Supply <system=...> variable in your environment' && tinker_abort=.true.

echo '*** Run Respa Integrator dynamic'
key_f=${system}_r
run_tinker 1 2500 2   # One   GPU  Execution
run_tinker 2 2500 2   # Two   GPUs Execution
run_tinker 4 2500 2   # Three GPUs Execution
run_tinker 8 2500 2   # Eight GPUs Execution
wait

key_f=${system}_rp
if [ -f ${key_f}.key ]; then
   run_tinker 2 2500 2   # Two   GPUs Execution
   run_tinker 4 2500 2   # Three GPUs Execution
   wait
fi

key_f=${system}_rp4
if [ -f ${key_f}.key ]; then
   run_tinker 8 2500 2   # Two   GPUs Execution
   wait
fi

echo '*** Run Respa1 Integrator dynamic'
key_f=${system}_R
run_tinker 1 2500 10   # One   GPU  Execution
run_tinker 2 2500 10   # Two   GPUs Execution
run_tinker 4 2500 10   # Three GPUs Execution
run_tinker 8 2500 10   # Eight GPUs Execution
wait

echo '*** Extract Simulation Speed'
[ -n "$system" ] && fetch_simulation_speed out_*.log > speed_${system}.txt

#save_logs bench_logs

echo !! $system system Bench Done !!

unset exe key_f tinker_abort
