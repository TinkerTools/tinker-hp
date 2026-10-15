!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module status_types  -- passive status declarations       ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
module status_types
   use iso_c_binding, only: c_int
   implicit none

   enum, bind(C)
      enumerator :: status_ok = 0
      enumerator :: status_missing_argument = 1
      enumerator :: status_invalid_argument = 2
      enumerator :: StatusDeviceError = 3
   end enum

   type :: status_t
      integer(c_int)     :: code    = status_ok
      character(len=240) :: message = ''
   end type status_t
end module status_types
