!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module status  --  explicit interoperable operation status  ##
!     ##                                                              ##
!     ##################################################################
#include "tinker_macro.h"

module status
   use iso_c_binding, only: c_int
   implicit none

   enum, bind(C)
      enumerator :: status_ok = 0
      enumerator :: status_missing_argument = 1
      enumerator :: status_invalid_argument = 2
   end enum

   type :: status_t
      integer(c_int)     :: code    = status_ok
      character(len=240) :: message = ''
   end type status_t

   interface
      module subroutine set_status(status_value, code, message)
         type(status_t), intent(out) :: status_value
         integer(c_int), intent(in)  :: code
         character(*)  , intent(in)  :: message
      end subroutine set_status

      module subroutine enrich_status(status_value, context)
         type(status_t), intent(inout) :: status_value
         character(*)  , intent(in)    :: context
      end subroutine enrich_status
   end interface
contains

   subroutine clear_status(status_value)
      implicit none
      type(status_t), intent(out) :: status_value

      status_value%code = status_ok
      status_value%message = ''
   end subroutine clear_status

   logical function failed(status_value)
      implicit none
      type(status_t), intent(in) :: status_value
      failed = status_value%code /= status_ok
   end function failed
end module status
