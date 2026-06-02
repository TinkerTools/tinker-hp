!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine command  --  get any command line arguments  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "command" uses the standard Unix-like iargc/getarg routines
!     to get the number and values of arguments specified on the
!     command line at program runtime
!
!
#include "tinker_precision.h"
subroutine command
   use argue
   implicit none
   integer i,iargc
   character*1 letter
   character*20 blank
!
!
!     initialize command line arguments as blank strings
!
   narg = 0
   blank = '                    '
   do i = 0, maxarg
      arg(i) = blank//blank//blank
   end do
!
!     get the number of arguments and store each in a string
!
   narg = command_argument_count ()
   if (narg .gt. maxarg)  narg = maxarg
   do i = 0, narg
      call get_command_argument (i,arg(i))
   end do
!
!     mark the command line options as unuseable for input
!
   listarg(0) = .false.
   do i = 1, narg
      listarg(i) = .true.
   end do
   do i = 1, narg
      letter = arg(i)(1:1)
      if (letter .eq. '-') then
         letter = arg(i)(2:2)
         call upcase (letter)
         if (letter.ge.'A' .and. letter.le.'Z') then
            listarg(i) = .false.
            listarg(i+1) = .false.
         end if
      end if
   end do
   return
end
