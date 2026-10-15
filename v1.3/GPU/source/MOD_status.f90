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
   use status_types, only: status_t, status_ok, status_missing_argument, status_invalid_argument, StatusDeviceError
   implicit none

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
