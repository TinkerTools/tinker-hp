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
!> @brief 
!> terminates execution due to a user request, a severe
!> error or some other nonstandard condition
!> @param no params
subroutine fatal
   use domdec
   use iounit
   implicit none
!
!
!     print a final warning message, then quit
!
   if (rank.eq.0) write (iout,10)
10 format (/,' TINKER is Unable to Continue; Terminating',&
   &' the Current Calculation',/)
   stop
end
