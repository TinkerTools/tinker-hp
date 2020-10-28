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
c     debug     logical flag to turn on full debug printing
c     holdup    logical flag to wait for carriage return on exit
c     abort     logical flag to stop execution at next chance
c
c     deb_Path   logical flag to print tinker Path
c     deb_Force  logical flag to print Forces information
c     deb_Energy logical flag to print Energy
c     deb_Polar  logical flag for polarisation potent
c
      module inform
      implicit none
      integer digits,iprint
      integer iwrite,isend
      logical silent,verbose
      logical debug,holdup,abort

      ! Debug static counter
      enum,bind(C)
        enumerator tindPath,tindForce,tindEnergy,tindAtom
        enumerator tindPolar,tindGdb
      end enum
      logical deb_Path,deb_Force,deb_Energy,deb_Atom,deb_Polar

      ! Inform separated Subroutines
      interface
        module subroutine initDebugEnv
        end subroutine
      end interface

      interface
        module subroutine info_minmax_pva
        end subroutine
      end interface
      interface
         module subroutine info_dyn()
         end subroutine
      end interface

      end
