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
      subroutine fatal
      use domdec
      use iounit
      implicit none
c
c
c     print a final warning message, then quit
c
      if (rank.eq.0) write (iout,10)
   10 format (/,' TINKER is Unable to Continue; Terminating',
     &           ' the Current Calculation',/)
      stop
      end
