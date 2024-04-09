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
c     dd_wverbose   logical flag to turn on extra information printing about dd
c     verbosepmetime   logical flag to turn on extra timings of PME printing
c     verboseforcestime   logical flag to turn on extra timings of forces
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
      logical silent,verbose,dd_verbose
      logical verbosepmetime,verboseforcestime
      logical debug,holdup,abort

      ! All program list
      enum,bind(C)
        enumerator analyze_a
        enumerator bar_a
        enumerator dynamic_a
        enumerator minimize_a
        enumerator testgrad_a
      end enum
      integer:: app_id=dynamic_a ! Only to be modifed inside a program

      ! Debug static counter
      enum,bind(C)
        enumerator tindPath,tindForce,tindEnergy,tindAtom
        enumerator tinMem,tindGdb,tindPolar
      end enum
      logical deb_Path,deb_Force,deb_Energy,deb_Atom,deb_Polar

c      integer dint1,dint2
c      integer,allocatable::dibuff(:)

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
c#ifdef USE_DETERMINISTIC_REDUCTION
c         module subroutine minmaxonef( vector,sz,name,mi_,ma_,on_ )
c         implicit none
c         integer sz
c         mdyn_rtyp vector(*)
c         character(*),optional,intent(in)::name
c         real(8)     ,optional,intent(inout):: mi_,ma_,on_
c!DIR$ ignore_tkr (r) vector
c         end subroutine
c#endif
         module subroutine minmaxonei( vector,sz,name )
         implicit none
         integer sz
         integer vector(*)
         character(*),optional,intent(in)::name
         end subroutine
         module subroutine minmaxonet( vector,sz,name )
         implicit none
         integer sz
         real*8 vector(*)
         character(*),optional,intent(in)::name
         end subroutine
c#if TINKER_MIXED_PREC
c         module subroutine minmaxoner( vector,sz,name )
c         implicit none
c         integer sz
c         real(r_p) vector(*)
c         character(*),optional,intent(in)::name
c!DIR$ ignore_tkr (r) vector
c         end subroutine
c#endif
      end interface

      interface normp
         module subroutine normt( array,n,val,p )
         integer siz
         integer,optional::p
         real*8 array(*)
         real*8 val
         end subroutine
c#if TINKER_MIXED_PREC
c         module subroutine normr( array,n,val,p )
c         integer siz
c         integer,optional::p
c         real(r_p) array(*)
c         real(r_p) val
c         end subroutine
c#endif
cc#if USE_DETERMINISTIC_REDUCTION
c         module subroutine normf( array,n,val,p )
c         integer siz
c         integer,optional::p
c         mdyn_rtyp array(*)
c         real(r_p) val
c!DIR$ ignore_tkr (r) array
c         end subroutine
c#endif
      end interface

c      interface
c      module subroutine set_dumpdyn_freq
c      end subroutine
c      end interface
c
c      interface
c      module subroutine check_loc(queue)
c      integer queue
c      end subroutine
c      end interface
c
c      interface
c      module subroutine check_loc1
c      end subroutine
c      end interface

      end
