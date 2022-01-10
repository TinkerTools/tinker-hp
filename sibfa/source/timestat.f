c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c   module timestat : time spent doing the different computations during a dynamic
c
c
      module timestat
      implicit none
      real*8 timeclear,timereneig,timecommstep,timeparam
      real*8 timeforcescomm,timedirreccomm,timereciptorquescomm
      real*8 timebondedvdw,timebonded,timenonbonded
      real*8 timenl
      real*8 timereal,timerec,timerecreccomm
      real*8 timerealdip,timerecdip
      real*8 timegrid1,timeffts,timescalar,timegrid2
      real*8 timestep
      save
      end
