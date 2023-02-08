c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine fatal  --  terminate the program abnormally  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "fatal" terminates execution due to a user request, a severe
c     error or some other nonstandard condition
c
c
#include "tinker_precision.h"
      subroutine fatal
      use domdec,only:rank
      use iounit,only:iout
      use mpi
      use mdstate
      implicit none
      integer errorcode,ierr
c
      if (track_mds) call mds_prt
c
c     print a final warning message, then quit
c
      if (rank.eq.0) write (0,10)
   10 format (/,' TINKER is Unable to Continue; Terminating',
     &          ' the Current Calculation',/)
      flush(iout)
!$acc wait
      call sleep(1)
      call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
      stop
      end
c
c
c
      subroutine fatal_acc
!$acc routine
      use domdec,only:rank
      implicit none

      if (rank.eq.0) then
         print*,'\n TINKER is Unable to Continue; Terminating',
     &           ' the Current Calculation \n'
         !print*,'-----Press Crtl-C to stop the calculation-----'
         stop
      end if
      end
c
c
c
      subroutine fatal1(msg,index)
      implicit none
      character(*),intent(in):: msg
      integer     ,intent(in):: index

 14   format(2X,/,'FATAL ERROR Detected in ',A,":l",I0,/)
      write(0,14) msg,index
      call fatal
      end subroutine
c
c
c
      subroutine fatal_device(msg)
      implicit none
      character(*),intent(in):: msg

 14   format(2X,'FATAL ERROR -- ',A,/,
     &       4x,'Unavailable feature for device platform')
      write(0,14) msg
      call fatal
      end
