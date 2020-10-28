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
      implicit none
      integer errorcode,ierr
c
c
c     print a final warning message, then quit
c
      if (rank.eq.0) write (iout,10)
   10 format (/,' TINKER is Unable to Continue; Terminating',
     &          ' the Current Calculation',/)
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
