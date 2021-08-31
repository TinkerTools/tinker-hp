c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c   module timestat : time spent doing the different computations during a dynamic
c     timestep: time spent for a complete time step
c     timeinte: time spent to update positions/velocities/accelerations
c     timereneig: time spent during "reneighboring" (deal with atoms
crossing domains)
c     timecommpos: time spent to communicate positions
c     timepara: time spent to regenerate local parameters
c     timegrad: time spent in gradient computation routines
c     timered: time spent in reduction operations to get energies
c     timetp: time spent in temperature/pressure control
c     timecommforces: time spent to communicate forces
c     timenl: time spent in neighbor list routines
c     timebonded: time spent in bonded forces routines
c     timevdw: time spent in vdw forces routines
c     timeelec: time spent in electrostatic forces routines
c     timepolar: time spent in polarization forces routines
c     timecleargrad: time spent to zero/sum forces arrays
c     
c     timereal: time spent to compute real space interactions (permanent)
c     timerec: time spent to compute reciprocal space interactions (permanent)
c     timecommforcesrec: time spent to communicate reciprocal forces
c     timecommforcesreal: time spent to communicate real space forces
c     timegrid: time spent to fill pme grid (permanent)
c     timefft: time spent to do ffts (permanent)
c     timescalar: time spent to do scalar product of pme (permanent)
c     timegrid2: time spent to get reciprocal space potential (permanent)
c     timerecreccomm: time spent in rec-rec communications
c 
c     
c
c
      module timestat
      implicit none
      real*8 timestep
      real*8 timeinte,timereneig
      real*8 timecommpos,timeparam
      real*8 timegrad,timered
      real*8 timetp,timecommforces
      real*8 timenl
      real*8 timebonded,timevdw,timeelec,timepolar
      real*8 timecleargrad
      real*8 timereal,timerec
      real*8 timecommforcesrec,timecommforcesreal
      real*8 timegrid,timefft,timescalar,timegrid2
      real*8 timerecreccomm
      save
      end
