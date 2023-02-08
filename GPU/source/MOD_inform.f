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
c     app_id    holds the running program identity
c     digits    decimal places output for energy and coordinates
c     iprint    steps between status printing (0=no printing)
c     iwrite    steps between coordinate dumps (0=no dumps)
c     isend     steps between socket communication (0=no sockets)
c     n_fwriten increments every time we force a frame to be written
c     monte_naccept  records the number of application of montecarlo barostat
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
c     tinEssai   variable for testing and developping
c
c     dint1 dint2 dibuff  useful for storing purpose

#include "tinker_macro.h"
      module inform
      implicit none
      integer digits,iprint
      integer iwrite,isend,n_fwriten
      logical silent,verbose
      logical debug,holdup,abort
      integer tinEssai
      integer mtc_nacc

      ! All program list
      enum,bind(C)
        enumerator analyze_a
        enumerator bar_a
        enumerator dynamic_a
        enumerator minimize_a
        enumerator radial_a
        enumerator testgrad_a
        enumerator pimd_a
      end enum
      integer:: app_id=dynamic_a ! Only to be modifed inside a program

      ! Debug static counter
      enum,bind(C)
        enumerator tindPath,tindForce,tindEnergy,tindAtom
        enumerator tinMem,tindGdb,tindPolar
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
        module subroutine info_minmax_pva(opt)
        integer,optional::opt
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
         module subroutine minmaxonei( vector,sz,name )
         implicit none
         integer sz
         integer vector(*)
         character(*),optional,intent(in)::name
!DIR$ ignore_tkr (r) vector
         end subroutine
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

      interface normp
         module subroutine normt( array,n,val,p )
         integer siz
         integer,optional::p
         real(t_p) array(*)
         real(r_p) val
!DIR$ ignore_tkr (r) array
         end subroutine
#if TINKER_MIXED_PREC
         module subroutine normr( array,n,val,p )
         integer siz
         integer,optional::p
         real(r_p) array(*)
         real(r_p) val
!DIR$ ignore_tkr (r) array
         end subroutine
#endif
#if USE_DETERMINISTIC_REDUCTION
         module subroutine normf( array,n,val,p )
         integer siz
         integer,optional::p
         mdyn_rtyp array(*)
         real(r_p) val
!DIR$ ignore_tkr (r) array
         end subroutine
#endif
      end interface

      interface
      module subroutine check_loc(queue)
      integer queue
      end subroutine
      end interface
      interface
      module subroutine check_loc1
      end subroutine
      end interface

      end
