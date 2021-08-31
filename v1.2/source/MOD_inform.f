c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module inform  --  control values for I/O and program flow  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     digits    decimal places output for energy and coordinates
c     iprint    steps between status printing (0=no printing)
c     iwrite    steps between coordinate dumps (0=no dumps)
c     isend     steps between socket communication (0=no sockets)
c     silent    logical flag to turn off all information printing
c     verbose   logical flag to turn on extra information printing
c     verbosepmetime   logical flag to turn on extra timings of PME printing
c     verboseforcestime   logical flag to turn on extra timings of forces
c     debug     logical flag to turn on full debug printing
c     holdup    logical flag to wait for carriage return on exit
c     abort     logical flag to stop execution at next chance
c
c
      module inform
      implicit none
      integer digits,iprint
      integer iwrite,isend
      logical silent,verbose
      logical verbosepmetime,verboseforcestime
      logical debug,holdup,abort
      save
      end
