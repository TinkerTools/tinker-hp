!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!
module timestat
   implicit none
   real*8 :: timestep !<time spent for a complete time step
   real*8 :: timeinte !<time spent to update positions/velocities/accelerations
   real*8 :: timereneig !<time spent during "reneighboring" (deal with atoms crossing domains)
   real*8 :: timecommpos !<time spent to communicate positions
   real*8 :: timeparam !<time spent to regenerate local parameters
   real*8 :: timegrad !<time spent in gradient computation routines
   real*8 :: timered !<time spent in reduction operations to get energies
   real*8 :: timetp !<time spent in temperature/pressure control
   real*8 :: timecommforces !<time spent to communicate forces
   real*8 :: timenl !<time spent in neighbor list routines
   real*8 :: timebonded !<time spent in bonded forces routines
   real*8 :: timevdw !<time spent in vdw forces routines
   real*8 :: timeelec !<time spent in electrostatic forces routines
   real*8 :: timepolar !<time spent in polarization forces routines
   real*8 :: timecleargrad !<time spent to zero/sum forces arrays
   real*8 :: timereal !<time spent to compute real space interactions (permanent)
   real*8 :: timerec !<time spent to compute reciprocal space interactions (permanent)
   real*8 :: timecommforcesrec !<time spent to communicate reciprocal forces
   real*8 :: timecommforcesreal !<time spent to communicate real space forces
   real*8 :: timegrid !<time spent to fill pme grid (permanent)
   real*8 :: timefft !<time spent to do ffts (permanent)
   real*8 :: timescalar !<time spent to do scalar product of pme (permanent)
   real*8 :: timegrid2 !<time spent to get reciprocal space potential (permanent)
   real*8 :: timerecreccomm !<time spent in rec-rec communications
   save
end
