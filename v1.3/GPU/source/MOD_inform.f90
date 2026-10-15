!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module inform  --  control values for I/O and program flow  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     app_id    holds the running program identity
!     digits    decimal places output for energy and coordinates
!     iprint    steps between status printing (0=no printing)
!     iwrite    steps between coordinate dumps (0=no dumps)
!     isend     steps between socket communication (0=no sockets)
!     idumpdyn  steps between restart dumps
!     n_fwriten increments every time we force a frame to be written
!     monte_naccept  records the number of application of montecarlo barostat
!     silent    logical flag to turn off all information printing
!     verbose   logical flag to turn on extra information printing
!     debug     logical flag to turn on full debug printing
!     holdup    logical flag to wait for carriage return on exit
!     abort     logical flag to stop execution at next chance
!
!     deb_Path   logical flag to print tinker Path
!     deb_Force  logical flag to print Forces information
!     deb_Energy logical flag to print Energy
!     deb_Polar  logical flag for polarisation potent
!     tinEssai   variable for testing and developping
!
!     dint1 dint2 dibuff  useful for storing purpose

#include "tinker_macro.h"
module inform
   use iso_c_binding , only: c_int64_t, c_int8_t
   implicit none
   integer digits,iprint
   integer iwrite,isend,n_fwriten,idumpdyn
   logical silent,verbose,dd_verbose
   logical debug,holdup,abort
   logical final_message_suppressed
   integer tinEssai
   integer mtc_nacc
   character(len=16) :: debug_phase = 'other'
   integer(c_int64_t):: debug_call  = 0_c_int64_t

   ! All program list
   enum,bind(C)
      enumerator analyze_a
      enumerator bar_a
      enumerator dynamic_a
      enumerator minimize_a
      enumerator radial_a
      enumerator testgrad_a
      enumerator pimd_a
      enumerator analyze_beads_a
      enumerator pibar_a
   end enum
   integer:: app_id=dynamic_a ! Only to be modifed inside a program

   ! Debug static counter
   enum,bind(C)
      enumerator tindPath,tindForce,tindEnergy,tindAtom
      enumerator tinMem,tindGdb,tindPolar
   end enum
   logical deb_Path,deb_Force,deb_Energy,deb_Atom,deb_Polar
   character*128 warning,sugest_vdw
!$acc declare create(warning,sugest_vdw)

   integer dint1,dint2
   integer,allocatable::dibuff(:)

   ! Inform separated Subroutines
   interface
      module subroutine initDebugEnv
      end subroutine
   end interface

   interface Tinker_shellEnv
      module subroutine Tinker_shellEnv_i(envName,envValue,defVal)
         character(*),intent(in):: envName
         integer,intent(in) :: defVal
         integer,intent(out):: envValue
      end subroutine
      module subroutine Tinker_shellEnv_f4(envName,envValue,defVal)
         character(*),intent(in):: envName
         real(4),intent(in) :: defVal
         real(4),intent(out):: envValue
      end subroutine
      module subroutine Tinker_shellEnv_f8(envName,envValue,defVal)
         character(*),intent(in):: envName
         real(8),intent(in) :: defVal
         real(8),intent(out):: envValue
      end subroutine
   end interface

   interface
      module subroutine set_warning
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
         integer n
         integer,optional::p
         real(t_p) array(*)
         real(r_p) val
!DIR$ ignore_tkr (r) array
      end subroutine
#if TINKER_MIXED_PREC
      module subroutine normr( array,n,val,p )
         integer n
         integer,optional::p
         real(r_p) array(*)
         real(r_p) val
!DIR$ ignore_tkr (r) array
      end subroutine
#endif
#if USE_DETERMINISTIC_REDUCTION
      module subroutine normf( array,n,val,p )
         integer n
         integer,optional::p
         mdyn_rtyp array(*)
         real(r_p) val
!DIR$ ignore_tkr (r) array
      end subroutine
#endif
   end interface

   interface
      module subroutine set_dumpdyn_freq
      end subroutine
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

   interface hash_me
      module procedure hash_real4
      module procedure hash_real8
      module procedure hash_int
      module procedure hash_int64
   end interface

   contains

   subroutine set_debug_phase(phase)
      character(*), intent(in) :: phase
      integer nchar
   
      nchar = min(len_trim(phase),len(debug_phase))
      debug_phase = '                '
      if (nchar > 0) then
         debug_phase(1:nchar) = phase(1:nchar)
      else
         debug_phase = 'other'
      endif
   end subroutine set_debug_phase

   subroutine hash_update(hash,bytes)
      integer(c_int64_t), intent(inout) :: hash
      integer(c_int8_t), intent(in) :: bytes(:)
      integer j
   
      do j = 1, size(bytes)
         hash = ieor(hash,iand(int(bytes(j),c_int64_t),255_c_int64_t))
         hash = hash * 1099511628211_c_int64_t
      end do
   end subroutine hash_update

   integer(c_int64_t) function hash_empty()
      hash_empty = -3750763034362895579_c_int64_t
   end function hash_empty
   
   integer(c_int64_t) function hash_real4(vector,count)
      real(4), intent(in) :: vector(*)
      integer, intent(in) :: count
      integer(c_int8_t) bytes(storage_size(vector(1))/8)
      integer(c_int64_t) hash
      integer i

      hash = hash_empty()
      do i = 1, count
         bytes = transfer(vector(i),bytes)
         call hash_update(hash,bytes)
      end do
      hash_real4 = hash
   end function hash_real4

   integer(c_int64_t) function hash_real8(vector,count)
      real(8), intent(in) :: vector(*)
      integer, intent(in) :: count
      integer(c_int8_t) bytes(storage_size(vector(1))/8)
      integer(c_int64_t) hash
      integer i

      hash = hash_empty()
      do i = 1, count
         bytes = transfer(vector(i),bytes)
         call hash_update(hash,bytes)
      end do
      hash_real8 = hash
   end function hash_real8

   integer(c_int64_t) function hash_int(vector,count)
      integer, intent(in) :: vector(*)
      integer, intent(in) :: count
      integer(c_int8_t) bytes(storage_size(vector(1))/8)
      integer(c_int64_t) hash
      integer i

      hash = hash_empty()
      do i = 1, count
         bytes = transfer(vector(i),bytes)
         call hash_update(hash,bytes)
      end do
      hash_int = hash
   end function hash_int

   integer(c_int64_t) function hash_int64(vector,count)
      integer(8), intent(in) :: vector(*)
      integer, intent(in) :: count
      integer(c_int8_t) bytes(storage_size(vector(1))/8)
      integer(c_int64_t) hash
      integer i

      hash = hash_empty()
      do i = 1, count
         bytes = transfer(vector(i),bytes)
         call hash_update(hash,bytes)
      end do
      hash_int64 = hash
   end function hash_int64

end
