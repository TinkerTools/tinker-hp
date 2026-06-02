!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine fatal  --  terminate the program abnormally  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "fatal" terminates execution due to a user request, a severe
!     error or some other nonstandard condition
!
!
#include "tinker_precision.h"
subroutine fatal
   use domdec,only:rank
   use iounit,only:iout
   use mpi
   use mdstate
   implicit none
   integer errorcode,ierr
!
   if (track_mds) call mds_prt
!
!     print a final warning message, then quit
!
   if (rank.eq.0) write (0,10)
10 format (/,' TINKER is Unable to Continue; Terminating',&
      &' the Current Calculation',/)
   flush(iout)
!$acc wait
   call sleep(1)
   call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
   stop
end
!
!
!
subroutine fatal_acc
!$acc routine
   use domdec,only:rank
   implicit none

   if (rank.eq.0) then
      print*,'\n TINKER is Unable to Continue; Terminating',&
         &' the Current Calculation \n'
      !print*,'-----Press Crtl-C to stop the calculation-----'
      stop
   end if
end
!
!
!
subroutine fatal1(msg,index)
   implicit none
   character(*),intent(in):: msg
   integer     ,intent(in):: index

14 format(2X,/,'FATAL ERROR Detected in ',A,":l",I0,/)
   write(0,14) msg,index
   call fatal
end subroutine
!
!
!
subroutine fatal_device(msg)
   implicit none
   character(*),intent(in):: msg

14 format(2X,'FATAL ERROR -- ',A,/,&
      &4x,'Unavailable feature for device platform')
   write(0,14) msg
   call fatal
end
