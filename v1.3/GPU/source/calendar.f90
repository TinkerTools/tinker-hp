!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine calendar  --  find the current date and time  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "calendar" returns the current time as a set of integer values
!     representing the year, month, day, hour, minute and second
!
!     note only one of the various implementations below should
!     be activated by removing comment characters
!
!
#include "tinker_precision.h"
subroutine calendar (year,month,day,hour,minute,second)
   implicit none
   integer year,month
   integer day,hour
   integer minute,second
!
!
!     use the standard "date_and_time" intrinsic function
!
   integer values(8)
   character*5 zone
   character*8 date
   character*10 time
   call date_and_time (date,time,zone,values)
   year = values(1)
   month = values(2)
   day = values(3)
   hour = values(5)
   minute = values(6)
   second = values(7)
!
!     use the obsolete "itime" and "idate" intrinsic functions
!
!     integer hms(3)
!     call itime (hms)
!     hour = hms(1)
!     minute = hms(2)
!     second = hms(3)
!     call idate (month,day,year)
   return
end
