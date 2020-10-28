run_tinker(){
   [ -z "$1" ] && echo 'Required Process MPI number' && return 1
   [ -z "$2" ] && echo 'Required number of iterations' && return 2
   [ -z "$3" ] && echo 'Required timestep' && return 2
   [ -z "$exe" ] && echo "Required executable's <exe=...> name to perform" && return 3
   [ -z "$key_f" ] && echo "Required keyfile <key_f=...> to perform" && return 3
   [ "$tinker_abort" = ".true." ] && return 4
#$1 corresponds to the number of mpi process in the dynamic
#$2 simulation timestep in [femtoseconds]
   local tinker_save=1000            # save frequency in picosecond. default to 1 ns of simulation
   local tinker_mode=2               # Simulation mode  [Keep it that way]
   local tinker_temperature=300      # Stabilisation temperature [K]
   set -x
   mpirun -np $1 $exe $system $2 $3 $tinker_save $tinker_mode $tinker_temperature -k $key_f > out_${key_f}_${1}gpu.log 2>out_${key_f}_${1}gpu.e
#Uncomment to enqueue tasks but do not forget to comment previous line
#   mpirun -np $1 $exe $system $2 $3 $tinker_save $tinker_mode $tinker_temperature -k $key_f > out_${key_f}_${1}gpu.log 2>out_${key_f}_${1}gpu.e &
   set +x
}

fetch_simulation_speed(){
   [ "$tinker_abort" = ".true." ] && return 4
   [ -z "$1" ] && echo 'Require filename to process' && return 1
   for i in $@ ; do
      echo
      echo "  >>>>>>  $i  <<<<<  "
      grep 'day' $i
   done
}

save_logs(){
   [ "$tinker_abort" = ".true." ] && return 4
   [ -z "$1" ] && echo "Required saving directory name" && return 1
   mkdir -p $1
   set -x
   mv out_* speed_* $1
   set +x
}

clean_logs(){
   [ "$tinker_abort" = ".true." ] && return 4
   set -x
   rm -fr out_* speed_*
   set +x
}
