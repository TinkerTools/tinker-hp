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
c     dint1 dint2 dibuff  useful for storing purpose

#include "tinker_precision.h"
#include "tinker_types.h"

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

      integer dint1,dint2
      integer,allocatable::dibuff(:)

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

      interface minmaxone
#ifdef USE_DETERMINISTIC_REDUCTION
         module subroutine minmaxonef( vector,sz,name,mi_,ma_,on_ )
         implicit none
         integer sz
         mdyn_rtyp vector(*)
         character(*),optional,intent(in)::name
         real(8)     ,optional,intent(inout):: mi_,ma_,on_
!DIR$ ignore_tkr (r) vector
         end subroutine
#endif
         module subroutine minmaxonet( vector,sz,name )
         implicit none
         integer sz
         real(t_p) vector(*)
         character(*),optional,intent(in)::name
!DIR$ ignore_tkr (r) vector
         end subroutine
#if TINKER_MIXED_PREC
         module subroutine minmaxoner( vector,sz,name )
         implicit none
         integer sz
         real(r_p) vector(*)
         character(*),optional,intent(in)::name
!DIR$ ignore_tkr (r) vector
         end subroutine
#endif
      end interface

      end
