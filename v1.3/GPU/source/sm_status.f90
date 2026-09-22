!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
#include "tinker_macro.h"

submodule (status) status_operations
contains

   module subroutine set_status(status_value, code, message)
      implicit none
      type(status_t), intent(out) :: status_value
      integer(c_int), intent(in) :: code
      character(*), intent(in) :: message

      status_value%code = code
      status_value%message = message
   end subroutine set_status

   module subroutine enrich_status(status_value, context)
      implicit none
      type(status_t), intent(inout) :: status_value
      character(*), intent(in) :: context

      if (len_trim(context) == 0) return
      if (len_trim(status_value%message) == 0) then
         status_value%message = context
      else
         status_value%message = trim(context)//': '//trim(status_value%message)
      end if
   end subroutine enrich_status

end submodule status_operations
