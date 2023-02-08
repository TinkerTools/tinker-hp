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
$> ci/install.sh nd
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
chk_gnuroot(){
if [[ -z ${GNUROOT+x} ]]; then
   GNUROOT=`which g++ | sed 's_/bin/g++__'`
   export GNUROOT
   cat << EOP
   !!! WARNING GNUROOT variable not find in environment !!!
   setting it to "$GNUROOT" according to g++ compiler location

   *** You might later need to export this variable for development purposes ***

   Installation will resume shortly
EOP
   sleep 4
fi
}
