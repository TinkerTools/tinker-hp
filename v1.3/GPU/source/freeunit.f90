!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ############################################################
!     ##                                                        ##
!     ##  function freeunit  --  gets an unopened logical unit  ##
!     ##                                                        ##
!     ############################################################
!
!
!     "freeunit" finds an unopened Fortran I/O unit and returns
!     its numerical value from 1 to 99; the units already assigned
!     to "input" and "iout" (usually 5 and 6) are skipped since
!     they have special meaning as the default I/O units
!
!
#include "tinker_precision.h"
function freeunit ()
   use iounit
   implicit none
   integer freeunit
   logical used
!
!
!     try each logical unit until an unopened one is found
!
   freeunit = 0
   used = .true.
   do while (used)
      freeunit = freeunit + 1
      if (freeunit.ne.input .and. freeunit.ne.iout) then
         if (freeunit .gt. 99) then
            write (iout,10)
10          format (/,' FREEUNIT  --  No Available Fortran',&
               &' I/O Units')
            call fatal
         end if
         inquire (unit=freeunit,opened=used)
      end if
   end do
   return
end
